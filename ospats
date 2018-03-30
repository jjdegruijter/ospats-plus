
# This file "ospats" contains the function ospats().
# Date: 2018-03-30
# Version: 1.01
# Author: Jaap de Gruijter

# CONTENT:
# 1. N x N matrix of generalized distances
# 2. Initial stratification
# 3. Contributions from strata to objective function (O)
# 4. Transfer grid points
# 5. Sample sizes
# 6. Stratified random sampling
# 7. Output of final results

# For process monitoring, uncomment println() lines

function ospats()
println("------------------------------------------------")
println("-------------------- START FUNCTION OSPATS ---")

##### SECTION 1. N x N MATRIX OF GENERALIZED DISTANCES

println("----- Calculating matrix of generalized distances --")
d2 = zeros(N,N)
for i = 1:(N-1)
  for j = (i+1):N
    d2[i,j] = ((z_pred[i] - z_pred[j])^2)/R2 + (s2[i] + s2[j]).*(1 - exp.(-3*(sqrt.((x[i] - x[j])^2 + (y[i] - y[j])^2)/range)))
  end
end
d2 = d2 + d2'

TOTd2 = sum(d2)/2
ObarH1 = sqrt(TOTd2)/N
stratcy = Array{UInt64,1}(N)
stratbest = Array{UInt64,1}(N)
Hbest = H_min

global stratbest, nhbest, Hbest, ObarFinal, nbest, nh_l, Nh

###########################################################
 for H = H_max : -1 : H_min        #  Start optimization of H

println("----------------------------Number of strata  :  ",H)

cbObj = Array{Float64,1}(H)
cbObj = zeros(Float64,H,1)
TotTransf = 0

##### SECTION 2. INITIAL STRATIFICATION

  missing = N - H*floor(UInt64, N/H)
  A = collect(1:H)
  B = vcat(A,A)
  repeat = N/H -2
  for rep = 1:repeat
    B = vcat(B,A)
  end
  fillup = collect(1:missing)
  B = vcat(B,fillup)
  strat0 = Array{UInt64,1}(N)
  v = collect(1:N)
  w = v[randperm(rng, N)]
  for i = 1:N
    strat0[w[i]] = B[i]
  end

##### SECTION 3. CONTRIBUTIONS FROM STRATA TO OBJECTIVE FUNCTION

Sd2 = zeros(1,H)
for strat = 1:H
  Sd2[strat] = 0
   for i = 1:(N-1)
    if strat0[i] == strat
       for j = (i+1):N
        if strat0[j] == strat
          Sd2[strat] = Sd2[strat] + d2[i,j]
        end
      end
    end
  end
end
Sd2Init = Sd2
cbObj = sqrt.(Sd2)
O = sum(cbObj)
ObarInit = O/N

##### SECTION 4. TRANSFER GRID POINTS

stratcy = Array{UInt64,2}(N,1)
stratcy = strat0
TotTransf = 0
TotCycle = 0
for cycle = 1:maxcycle
  transfers = 0
  u = randperm(N)
  for t = u
    Delta = 0
    change = 0
    A = stratcy[t]
    ij = find(stratcy .== A)
    dA = sum(d2[t,ij])
    sumd2tinA = dA
    Sd2Amint = Sd2[A] - sumd2tinA
    cbObjA = sqrt.(abs.(Sd2Amint))
    for stratnr = 1:H
      Delta = 0
      sumd2plus = 0
      if stratnr != A
        B = stratnr
        ij = find(stratcy .== B)
        dB = sum(d2[t,ij])
        sumd2plus = dB
        cbObjB = sqrt.(abs.(Sd2[B] + sumd2plus))
        Delta = cbObjA + cbObjB - cbObj[A] -cbObj[B]
        if Delta < O*1e-10
          change = 1
          transfers = transfers + 1
          stratcy[t] = B            # update stratification
          Sd2[A] = Sd2[A] - sumd2tinA
          Sd2[B] = Sd2[B] + sumd2plus
          cbObj = sqrt.(abs.(Sd2))
          O = sum(cbObj)
          Delta = 0
        end                       # if Delta < Obj*1e-10
      end                         # if stratnr != A
      if change ==1 break end
    end                           # for strat=1:H
  end                             # for t=u

  # println("cycle ", cycle, "     transfers = ", transfers)
  TotTransf = TotTransf + transfers
  if transfers == 0 break end     # stopping rule
  TotCycle = cycle
end                               # for cycle=1:maxcycle

println("Total number of transfers = ", TotTransf)
println("Number of iteration cycles = ", TotCycle)
O = sum(cbObj)
ObarFinal = O/N

##### SECTION 5. SAMPLE SIZES

##### Subsection 5.1 Size of the grid-strata
Nh = Array{UInt64}(1,H)
Nh = zeros(H)
for h = 1:H
  k = find(stratcy .== h)
  Nh[h] = length(k)
end

##### Subsection 5.2 Total sample size before correction for
##### roundoff error.
n_pred = (CP*Area*Z_gamma*ObarFinal/(f*sqrt.(2)))^(2/3)

##### Subsection 5.3 Neyman allocation
sum_ahOh = 0
for h = 1:H
   ahOh = Nh[h]*cbObj[h]
   sum_ahOh = sum_ahOh + ahOh
end
nh = Array{Float64}(H)
nhbest = Array{Float64}(H)
for h = 1:H
  nh[h] = n_pred *Nh[h] *cbObj[h]/sum_ahOh
end
nh = round.(nh)

##### Subsection 5.4 Correction of total sample size to
##### avoid difference with sum of sizes in strata.
n_pred = sum(nh)
n_pred = convert(UInt64,n_pred)

##### Subsection 5.5 Check on smallest sample size in strata.
nh_low = indmin(nh)
nh_l = nh[nh_low]
nh_min = Array{Float64,1}
nh_min = nh_minim

##### Subsection 5.6 Update and output of intermediate results.
Hbest = H
stratbest = stratcy
nhbest = nh
nbest = n_pred

println("intermediate results :")
println("Sample size :   ", nbest)
println("Sample sizes in strata : ", nhbest)
println("Smallest sample size allocated to a stratum :  ", nh_l)

if nh_l >= nh_min break end      # stop lowering H

end                                         # end optimization of H

###################################################################

##### SECTION 6. STRATIFIED RANDOM SAMPLING WITH FINAL
##### STRATIFICATION AND NEYMAN ALLOCATION

 n_tot = nbest   # sum(nhbest)
 n_tot = round(n_tot)
 n_tot = convert(UInt64,n_tot)
 points = Array{UInt64}(1,n_tot)
 xs = Array{Float64}(1,n_tot)
 ys = Array{Float64}(1,n_tot)

 for h = 1:Hbest
  k = find(stratbest .== h)      # numbers of points in stratum h
  stratsize = length(k)
  v = collect(1:stratsize)
  w = v[randperm(rng,stratsize)]  # randomized indexes to points

  f=0
  for i = 1:nhbest[h]
    f=f+1                     # making Neyman allocations integer
  end

  k_rand = k[w]               # randomized points in stratum h
  points_h = k_rand[1:f]      # put the first f points in sample

  if h == 1                   # concatenate all H vectors of points
    points = points_h
  elseif h > 1
    points = vcat(points,points_h)
  end
 end
xs = x[points]                # get x coordinates of sample points
ys = y[points]                # get y coordinates of sample points
strata = stratbest[points]    # get stratum numbers of sample points
sampnr = collect(1:n_tot)     # make sample numbers

strs = DataFrame()
strs[:SampleNr] = sampnr
strs[:StratNr] = strata
strs[:PointNr] = points
strs[:X] = xs
strs[:Y] = ys

stratif = DataFrame()
stratif[:X] = x
stratif[:Y] = y
stratif[:Stratum] = stratbest

##### SECTION 8. OUTPUT OF FINAL RESULTS

CSV.write("stratification", stratif)
CSV.write("Sample", strs)
println("    ")
println("FINAL RESULTS: ")
println("Obar (O/N): ", ObarFinal)
println("Number of strata : ", Hbest)
println("Total sample size : ",nbest)
println("Sample sizes in strata : ", nhbest)
println("Smallest sample size allocated to a stratum :  ",  nh_l)
println("Sizes of grid-strata : ", Nh)

println("--------------------- END FUNCTION OSPATS ---")
end                         # function ospats
