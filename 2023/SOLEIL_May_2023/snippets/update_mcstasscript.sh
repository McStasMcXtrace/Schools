python3 -m venv my_environment
source my_environment/bin/activate
pip3 install --proxy=http://195.221.10.6:8080 mcstasscript jupyterlab ipympl
sed -i s+/Applications/McXtrace-1.5.app/Contents/Resources/mcxtrace/1.5/+/usr/share/mcxtrace/3.2-20230508+g my_environment/lib/python3.9/site-packages/mcstasscript/configuration.yaml
ipython3 kernel install --user --name=my_environment
