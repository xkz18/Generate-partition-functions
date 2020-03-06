def Partition(q_1,Nmax):
    
    Nmax = Nmax+1
    deriv_f1=np.zeros(Nmax)
    deriv_f2=np.zeros(Nmax)
    f3=np.zeros(Nmax)
    f2=np.zeros(Nmax)
    f2_new=np.zeros(Nmax)
    Qn=np.zeros(Nmax)
    fac=np.zeros(Nmax)
    ntot = Nmax-2

    fac[0]=1.0
    fac[1]=1.0
    for i in range(2,ntot+1):
        fac[i]=fac[i-1]*i
    
    for i in range(0,ntot+1):
        deriv_f1[i]= q_1[i+1]*(i+1)
    
    Qn[0]=1.0/fac[0]   
    f2=deriv_f1
    Qn[1]=f2[0]/fac[1]
    

    for i in range(2,ntot+1):
        deriv_f2=np.zeros(Nmax)
        f3=np.zeros(Nmax)
        f2_new=np.zeros(Nmax)


        for j in range(0,ntot+1):
            deriv_f2[j]=f2[j+1]*(j+1)

        for j in range(0,ntot+1):
            for k in range(0,ntot+1-j):
                m=j+k
                f3[m]=f3[m]+f2[j]*deriv_f1[k]
      
        for j in range(0,ntot+1):
            f2_new[j]=f3[j]+deriv_f2[j]
        
        f2=f2_new
        Qn[i]=f2[0]/fac[i]
   
    InvQn=1.0/Qn[ntot]
    #InvQn=1.0/Qn[ntot-1]
    n_output=np.zeros(Nmax)

    for i in range(0,ntot+1):
        n_output[i] = q_1[i]*InvQn*Qn[ntot-i]
    #print("InvQn",InvQn)
    #print("noutput",n_output)
    return InvQn,n_output