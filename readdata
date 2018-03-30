# This file "readdata" contains the function readdata().
# Date: 2018-03-30
# Version: 1.01
# Author: Jaap de Gruijter

function readdata()

global x, y, z_pred, s2, N
##################################### Reading the data file
# Data = readtable(filename, header=false, separator=',')
Data = CSV.read(filename; delim=',')
x = Data[:1]
y = Data[:2]
# nr = Data[:3]        grid point identifier
z_pred = Data[:4]
s2 = Data[:5]
N = length(x)
println("Grid size : ", N)

##################################### Checking process parameters

if maxcycle < 0                              # check maxcycle
  println("ERROR: maxcycle is negative")
  quit()
end
type_maxcycle = typeof(maxcycle)
if type_maxcycle != Int64
  print("ERROR: maxcycle is not an integer")
  quit()
end

if in < 1                                    # check in
  println("ERROR: in is smaller than 1")
  quit()
end
type_in = typeof(in)
if type_in != Int64
  print("ERROR: in is not an integer")
  quit()
end

if R2 < 0                                    # check R2
  println("ERROR: R2 is negative")
  quit()
end
if R2 > 1
  println("ERROR: R2 is larger than 1")
  quit()
end

if range <= 0                                 # check range
  println("ERROR: range is zero or negative")
  quit()
end

println("Seed for random number generator : ",seed)
println("Sampling interval :  ", in)
println("R2 :  ", R2)
println("Range :  ", range)
println("Maximum number of iteration cycles :  ", maxcycle)
println("Minimum sample size allowed in the strata : ", nh_minim)
println("Minimum number of strata : ", H_min)
println("Maximum number of strata : ", H_max)
println("Carbon offset price (Aus dollar per Mg) : ", CP)
println("Cost of obtaining data per grid point (Aus dollar) : ", f)
println("Surface area of the farm (ha) : ", Area)

end
