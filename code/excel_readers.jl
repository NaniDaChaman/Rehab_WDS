using ExcelReaders

data = readxl("input/Sample_input_jt.xls", "Node Data!A2:E10")
for i = 1:9 
    println(data[i,1])
end 
#print(data[:,1])
#data is read in as a 2d array !