//Sets input and output
input = "/Users/kristal/desktop/seq_smFISH/CROP_WT/";
output = "/Users/kristal/desktop/seq_smFISH/CROP_WT/crop/"

//Loop
list = getFileList(input);
for (i=0; i < list.length; i++) {
	if (startsWith(list[i], "hyb4"))
	if (endsWith(list[i], "39_MIP.tif"))
	action(input, output, list[i]);

function action(input, output, filename) {
	
open(input + filename); {

t=File.nameWithoutExtension;	


// Crop and save

setTool("polygon");
makePolygon(1614,17,1505,19,1403,51,1383,108,1395,201,1427,281,1483,419,1527,494,1552,457,1572,298,1657,191,1718,32);
run("Duplicate...", "duplicate");
saveAs("Tiff", output+t+"_6.tif");
close();

makePolygon(1653,1498,1668,1683,1960,1665,2001,1481,1908,1454,1695,1461);
run("Duplicate...", "duplicate");
saveAs("Tiff", output+t+"_7.tif");
close();

makePolygon(1650,1552,1677,1726,1618,1845,1533,1916,1411,1920,1363,1880,1383,1765,1435,1653,1550,1572);
run("Duplicate...", "duplicate");
saveAs("Tiff", output+t+"_8.tif");
close();

//makePolygon(1370,1577,1465,1577,1550,1573,1585,1573,1568,1700,1479,1812,1404,1777,1347,1693,1357,1607);
//run("Duplicate...", "duplicate");
//saveAs("Tiff", output+t+"_5.tif");
//close();



close();
}
}
}
