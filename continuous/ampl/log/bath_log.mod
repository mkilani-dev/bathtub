
# MODEL PARAMETERS
param P    ;
param K    ;
param vf   ;
param M0   ;
param Mbar ;
param TMAX ;
param L    ;
param T0   := M0 / vf ; 

param n ;
param q ;

# APPROXIMATION GRID
param B {i in 1..n+q};
param q0 ;
# LOG UTILITY FUNCTION
# --> LOG
param R0 ;
param R1 ;
param S0 ;
param S1 ;

# A SMALL VALUE
param eps  ;


# GENERATE INITIAL SOLUTION WITH UNIFORM ENTRIES
#var e {i in 1..n} >= 0.0, :=
#      (B[i+1]-B[i])/(B[n+1]-B[1]) ; 
var e {1..n} >=0 ;

var k  { i in 1..n+q-1 }  # traffic density
       = sum {j in max(1,i-q+1)..min(i,n)} e[j];
       
var v  { i in 1..n+q-1 } # traffic speed
       = vf * (1 - P * k[i] / K );

var Dt { i in 1..n+q-1 } # travel time from B[i] to B[i+1]
       = ( B[i+1] - B[i] ) / v[i];

var t  { i in 1..n+q } = # clock time at B[i]
       if i = 1 then T0
       else t[i-1] + Dt[i-1] ;


# LOGARITHMIC UTILITY FUNCTION
var ue_log {i in 1..n}
       =  R0 * log( R1 * t[i] ) ;

var ux_log {i in 1..n}
       =  S0 * log( S1 * ( TMAX - t[i+q] ) ) ;

maximize util_log: sum{i in 1..n} P * ( ue_log[i] * e[i] + ux_log[i] * e[i] ) ;



# MINIMIZE THE SHORTEST END TIME OF THE RUSH HOUR
minimize last_arrival_time: t[n+q] ; 



# CONSTRAINTS OF THE PROBLEM 
  # 1. All agents enter the bathtub
subject to entry : sum{i in 1..n}  e[i] == 1.0 ;

  # 2. Jam capacity should not be reached
subject to capacity { i in q .. n } : k[i] <= K / P - eps ;


