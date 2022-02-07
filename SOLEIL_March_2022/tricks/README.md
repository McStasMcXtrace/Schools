## Fix the html doc installation ##

🚳️ The `mxdoc` command is not functional with McXtrace 1.5 installation, which prevents the MxGUI/Help menu items to work. To fix it:
<pre>
sudo apt-get install mcxtrace-tools-python-mccodelib-1.5 mcxtrace-tools-python-mxdoc-1.5
</pre>
followed by
<pre>
sudo mxdoc -i
</pre>

## Fix for the missing gfortran library

🚭️ On some systems, the pre-compiled version of `cif2hkl` given with McXtrace (and McStas) gives the error:
```
cif2hkl: error while loading shared libraries: libgfortran.so.5: cannot open shared object file: No such file or directory
```

To get rid of this, one just needs to recompile `cif2hkl`, e.g.:
``` bash
sudo apt install gfortran
cd /usr/share/mcxtrace/1.5/libs/cif2hkl
sudo gfortran -O2 -o cif2hkl cif2hkl.F90 -lm
sudo cp cif2hkl /usr/bin/
sudo rm *.mod
```

## Configuring the PGPLOT plotters for not using Giza
🚳🚭️ Ubuntu 18.04 enables the Giza backend when using our PGPLOT
mxplot/mxdisplay backends. To configure for the slightly more stable
legacy PGPLOT libraries, please

``` bash
cd /usr/lib/x86_64-linux-gnu
sudo ln -sf ../libpgplot.so libpgplot.so.0.0.0
sudo ln -sf ../libcpgplot.so libcpgplot.so.0.0.0
```

## Overview of mxplot and mxdisplay plotters
See [our wiki](https://github.com/McStasMcXtrace/McCode/wiki) and specifically the [mxplot](https://github.com/McStasMcXtrace/McCode/wiki/mcplot-variants---table-overview) and [mxdisplay](https://github.com/McStasMcXtrace/McCode/wiki/mcdisplay-variants---table-overview) pages
