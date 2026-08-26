def writeKoDotIn(s1, Y):
    import sys

    print('writeKoDotIn is running ...') 
    s1 = f"{s1:.3f}"
    Y = f"{Y:.1e}"

    if 'e+' in Y:
        base, exp = Y.split('e+')
        Y = f"{base}e+{int(exp)}"
    elif 'e-' in Y:
        base, exp = Y.split('e-')
        Y = f"{base}e-{int(exp)}"

    with open("ko.in", "w", newline="") as f:
        f.write('jGeometry:                          Mbar  cm/us   g/cc          g/cc  cm/us unitless            Megabar Megbar Megbar Mbar-cc/K-g' + '\n')                             
        f.write('iEOS iStr Nodes   Length   xstart    P_0    U_0  Rho_0    E_0 Rho_00      C     s1     s2  gamma Yeld,Y   mu,G  pfrac     Cv    a_p    a_e    p_e    p_s  rho_0 beta    n    b       AJO       BJO       CJO       MJO       NJO       TJO     Pos' + '\n')               
        f.write('   1    1   306    0.306   -0.306    0.0  0.039  16.69  0.000  0.000 0.3414  ' +  str(s1) +'   0.00   1.67 '+str(Y)+'   0.69  0.000 1.6e-5   0.00   0.00   0.00   0.00   0.00  0.0  0.0  0.0       0.0       0.0       0.0       0.0       0.0       0.0     0.3                             '+ '\n')                       
        f.write('   1    1   960    0.960   -0.000    0.0   0.00  16.69  0.000  0.000 0.3414  ' +  str(s1) +'   0.00   1.67 '+str(Y)+'   0.69  0.000 1.6e-5   0.00   0.00   0.00   0.00   0.00  0.0  0.0  0.0       0.0       0.0       0.0       0.0       0.0       0.0     0.3                             '+ '\n')                       
        f.write("""
Boundary Conditions:                                                                                                                                             
   iend    P_0    U_0  Rho_0    E_0    V_0                                                                                                                        
    -1 0.0e-1     .0   2.70  0.000    1.0                                                                                                                                  
     1 0.0e-1     .0   2.70  0.000    1.0                                                                                                                                      

tstop  dtskip  (micro seconds):
5.0    0.01                                                                                                                                                                        
Do not erase this statement
    1    1  1000    1.000   -1.000  .0E-6 0.0451  2.340  0.000  0.000 0.0551  4.520   0.00  2.000 0.0e-2  0.248  0.000 1.0e-9   0.00   0.00   0.00   0.00   0.00  0.0  0.0  0.0       0.0       0.0       0.0       0.0       0.0       0.0     0.3                       
    1    1   100    0.100    0.000  .0E-6 0.0000  8.930  0.000  0.000  0.390  1.490   0.00  1.990 0.0e-3  0.477 3.9e-5 1.0e-9   0.00   0.00   0.00   0.00   0.00  0.0  0.0  0.0       0.0       0.0       0.0       0.0       0.0       0.3     0.3
    """ + '\n')

if __name__ == "__main__":
    import sys

    # manually set values if running directly
    #s1 = 0.123
    #Y = 1.2e-7

    # Get command-line arguments
    if 's1' not in locals():
        s1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.2
    if 'Y' not in locals():
        Y = float(sys.argv[2]) if len(sys.argv) > 2 else 7.7e-3

    writeKoDotIn(s1, Y)