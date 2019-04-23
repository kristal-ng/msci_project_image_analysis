//Sets input and output
input = "/.../";
output = "/.../";

//Loop
list = getFileList(input);
for (i=0; i < list.length; i++) {
	if (endsWith(list[i], ".tif"))
	//if (startsWith(list[i], "hyb0"))
	action(input, output, list[i]);

function action(input, output, filename) {
	
open(input + filename); {

//Set grayscale
run("Channels Tool...");
Stack.setDisplayMode("grayscale");

//Split channels into channel 1, 2, 3, 4, 5
run("Split Channels");
t=File.nameWithoutExtension;	



//Saves each grayscale image for each channel with original name and channel number
selectWindow("C5-"+t+".tif");
saveAs("Tiff", output+t+"_cy5.tif");
close();
selectWindow("C4-"+t+".tif");
saveAs("Tiff", output+t+"_cy3.tif");
close();
selectWindow("C3-"+t+".tif");
close();
selectWindow("C2-"+t+".tif");
saveAs("Tiff", output+t+"_fid.tif");
close();
selectWindow("C1-"+t+".tif");
saveAs("Tiff", output+t+"_dapi.tif");
close();

}
}
}
