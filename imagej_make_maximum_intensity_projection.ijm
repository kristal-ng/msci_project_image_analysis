//Sets input and output
input = "~/input folder/";
output = "~/output folder/";

//Loop through the folder, open only .tif files  
list = getFileList(input);
for (i=0; i < list.length; i++) {
	if (endsWith(list[i], ".tif"))
	action(input, output, list[i]);

function action(input, output, filename) {
	
open(input + filename); {

//Make composite image from z-stack with maximum intensity projection
//Will trim 10 sections from the beginning and at the end

title = getTitle();
run("Brightness/Contrast...");
run("Slice Remover", "first=206 last=215 increment=1");
run("Slice Remover", "first=1 last=10 increment=1");
run("Deinterleave", "how=5 keep");
run("Merge Channels...", "c1=[" + title + " #2] c2=[" + title + " #3] c3=[" + title + " #1] c4=[" + title + " #4] c5=[" + title + " #5] create");
selectWindow("Composite");
Stack.setChannel(1);
setMinAndMax(1,5000);
Stack.setChannel(2);
setMinAndMax(1, 5000);
Stack.setChannel(3);
setMinAndMax(1,5000);
Stack.setChannel(4);
setMinAndMax(1,5000);
Stack.setChannel(5);
setMinAndMax(1,5000);

run("Channels Tool...");
run("Brightness/Contrast...");
run("Z Project...", "projection=[Max Intensity]");

//Saves max intensity projection window as .tif in output folder with "_MIP" tag after title
t=File.nameWithoutExtension;	
selectWindow("MAX_Composite");
saveAs("Tiff", output+t+"_MIP.tif");

close();
close();
close();
}
}
}
