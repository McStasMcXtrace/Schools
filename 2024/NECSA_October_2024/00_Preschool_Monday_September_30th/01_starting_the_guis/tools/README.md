# Tool correspondance, McStas vs McXtrace
<hr>
<table>
<tr><td><strong>McStas</strong></td><td><strong>McXtrace</strong></td></tr>
<tr><td>main UI<br><pre>mcgui</pre></td><td>main UI<br><pre>mxgui</pre></td></tr>
<tr><td>documentation<br><pre>mcdoc</pre></td><td>documentation<br><pre>mxdoc</pre></td></tr>
<tr><td>plot output data<br><pre>mcplot(-pyqtgraph)<br>mcplot-matplotlib<br>mcplot-html</pre></td>
	<td>plot output data<br><pre>mxplot(-pyqtgraph)<br>mxplot-matplotlib<br>mxplot-html</pre></td></tr>
<tr><td>visualize instrument<br><pre>mcdisplay(-pyqtgraph)<br>mcdisplay-webgl(-classic)<br>mcdisplay-matplotlib<br>mcdisplay-cad<br>mcdisplay-mantid</pre></td>
	<td>visualize instrument<br><pre>mxdisplay(-pyqtgraph)<br>mxdisplay-webgl(-classic)<br>mxdisplay-matplotlib<br>mxdisplay-cad<br>-</pre></td></tr>
<tr><td>"System" directory<br><pre>$MCSTAS<br>(%MCSTAS% on Windows)</pre></td><td>"System" directory<br><pre>$MCXTRACE<br>(%MCXTRACE% on Windows)</pre></td></tr>
<tr><td>Data file directory<br><pre>$MCSTAS/data<br>(%MCSTAS%\data on Windows)</pre></td><td>Data file directory<br><pre>$MCXTRACE/data<br>(%MCXTRACE%\data on Windows)</pre></td></tr>
<tr><td>example instrument directory<br><pre>$MCSTAS/examples<br>(%MCSTAS%\examples on Windows)</pre></td>
	<td>example instrument directory<br><pre>$MCXTRACE/examples<br>(%MCXTRACE%\examples on Windows)</pre></td></tr>
<tr><td>component categories:<br><pre>
sources
monitors
optics
samples
union
sasmodels
misc
contrib
obsolete
<br></pre></td><td>component categories:<br><pre>
sources
monitors
optics
samples
union
sasmodels
misc
astrox
contrib
obsolete
<br></pre></td></tr>
<tr><td>"share" directory<br><pre>$MCSTAS/share<br>(%MCSTAS%\share on Windows)</pre></td>
	<td>"share" directory<br><pre>$MCXTRACE/share<br>(%MCXTRACE%\share on Windows)</pre></td></tr>
</table>
<hr>