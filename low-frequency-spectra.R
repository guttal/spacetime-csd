# Written by Sabiha Majumder

# this code calls the functions powerspectrum and coarse_grain. 
#These function codes should be saved in the same folder

a=read.table("filename2.txt",header=FALSE)
a=a[,]
X=181;  # Number of driver-values in the transect 
n=128;
b=matrix(-1,n,n);   # Makes the matrix from the data
mat_cg=matrix(n/4,n/4);  # Coarse-grained matrix   
spectrum<-numeric(X);  # an array which saves the dft values

for (x in 1:X){
for ( i in 1:n){
  for (j in 1:n){
    b[i,j]=a[(x-1)*n*n+(i-1)*n+j];
  }
}

mat_cg=coarse_grain(b,4)
c=powerspectrum(b);
c1=mean(c[1:8])

spectrum[x]<-c1
}

