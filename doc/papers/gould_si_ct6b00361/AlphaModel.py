import numpy as np

"""
Python class for dealing with polarizability model data
described in [REF].

Methods are:
 
A=AlphaModel(ModelFileName) # Initialise model

[a1,O1,a2,O2]=A.GetModel((Z,N)) # Get model parameters as array [a1,O1,a2,O2]
Alpha(omega)=A.GetAlpha((Z,N), omega=omega) # Get alpha at omega

C6=A.GetC6( (Z,N) ) # Get same-species C6 coefficient
C6=A.GetC6( [(Z1,N1),(Z2,N2)] )  # Get different-species C6 coefficient
print "%8s"%( GetC6Str(self, (Z,N)) )  # Print formatted


"""

# For quick evaluation of C6 coefficients
def FC6(a1,o1,a2,o2):
    return 1.5*a1*a2/(o1+o2)/o1/o2

# For quick testing of pairs
def Pair(Z):
    if isinstance(Z[0],(list,tuple)) and (len(Z)>1):
        return Z[0],Z[1]
    else:
        return Z,Z

# Quick formatting
def QFmt(C6):
    if (C6<1.0):
        return "%8.3f"%(C6)
    elif (C6<10.0):
        return "%8.2f"%(C6)
    elif (C6<100.0):
        return "%8.1f"%(C6)
    else:
        return "%8.0f"%(C6)

"""
Class implementing the model
"""
class AlphaModel:

    """Initialisation.
    By default loads the benchmark set (must be present)"""
    def __init__(self, file="ModelPGG_Scaled.dat"):
        self.DataFile=file
        self.Data=np.loadtxt(file)
        self.Z=self.Data[:,0]
        self.N=self.Data[:,1]
        self.Model=self.Data[:,3:7]
        
        # Build an index
        self.Indx={}
        II=self.Z*1000 + self.N # There are < 1000 atoms and ions
        for k in range(len(II)):
            self.Indx[II[k]]=k

    """Get the model data for ZN=(Z,N) using the index"""
    def GetModel(self, ZN):
        II=ZN[0]*1000+ZN[1]
        if II in self.Indx:
            return self.Model[self.Indx[II]]
        else:
            raise Exception("Element not defined")

    """Use the model to calculate alpha at omega.
    omega can be an array."""
    def GetAlpha(self, ZN, omega=0.0):
        M=self.GetModel(ZN)
        return M[0]/(M[1]**2 + omega**2) + M[2]/(M[3]**2 + omega**2) 

    """Get the C6 coefficient using the model.
    Here ZN=(Z,N) returns a same-species coefficient
    and  ZN=( (Z1,N1), (Z2,N2) ) returns the cross-coefficient."""
    def GetC6(self, ZN):
        ZN1,ZN2=Pair(ZN)

        M1=self.GetModel(ZN1)
        M2=self.GetModel(ZN2)

        C6 =  FC6(M1[0], M1[1], M2[0], M2[1]) \
            + FC6(M1[0], M1[1], M2[2], M2[3]) \
            + FC6(M1[2], M1[3], M2[0], M2[1]) \
            + FC6(M1[2], M1[3], M2[2], M2[3])

        return C6

    """Same as GetC6 but returns formatted string of 8 characters."""
    def GetC6Str(self, ZN): return QFmt(self.GetC6(ZN))


if __name__ == "__main__":
    print "A=AlphaModel()"
    A=AlphaModel()
    print "A.GetModel((4,4))          =",(A.GetModel((4,4)))
    print "A.GetAlpha((4,4),omega=0.) =",(A.GetAlpha((4,4), omega=0.0))
    print "A.GetC6((4,4))             =",(A.GetC6((4,4)))
    print "A.GetC6(((2,2),(11,11)))   =",(A.GetC6(((2,2),(11,11))))
    print "A.GetC6Str(((2,2),(11,11)))=",(A.GetC6Str(((2,2),(11,11))))

