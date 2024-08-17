LIB("rinvar.lib"); 
LIB("invar.lib"); 
LIB("matrix.lib");

int N = 3;; 


ring R = 0, (g(1..N*N), a(1..2*N*N)), dp; 

matrix g[N][N]; 
matrix A1[N][N]; 
matrix A2[N][N]; 

matrix g_adjoint[N][N]; 

for(int i=1;i<=N;i++)
{
	for(int j=1;j<=N;j++)
	{
		g[i,j] = g((i-1)*N+j); 
		A1[i,j] = a((i-1)*N+j); 
		A2[i,j] = a(N+(i-1)*N+j); 
	}
	kill j; 
}
kill i; 

intvec r=1..N; 
for(int i=1;i<=N;i++)
{
	for(int j=1;j<=N;j++)
	{
		g_adjoint[i,j] = (-1)^(i+j)*det(submat(g,delete(r,j),delete(r,i))); 
	}
	kill j;
}
kill i;
kill r;

matrix a1 = g*A1*g_adjoint; matrix a2 = g*A2*g_adjoint; 

matrix m[2*N*N][2*N*N]; 



for(int i=1;i<=2*N*N;i++)
{
	for(int j=1;j<=2*N*N;j++)
	{
		poly r_temp, c_temp; 

		if(i<=N*N)
		{
			r_temp = a1[1+((i-1) div N),1+((i-1)%N)]; 
		}
		else
		{
			r_temp = a2[1+((i-N*N-1) div N),1+((i-N*N-1)%N)]; 
		}

		if(j<=N*N)
		{
			c_temp = `string(A1[1+((j-1) div N),1+((j-1)%N)])`; 
		}
		else
		{
			c_temp = `string(A2[1+((j-N*N-1) div N),1+((j-N*N-1)%N)])`; 
		}

		m[i,j] = diff(r_temp, c_temp); 
		kill r_temp;
		kill c_temp;
	}
	kill j;
}
kill i;

SL(N);
setring(Invar::group); 
map f = R, g(1..N*N); 
matrix m = f(m); 

// invar(m); 


/*
matrix m[8][8]; 
m[1,1]=g(1)*g(4); 
m[1,2]=-g(1)*g(3); 
m[1,3]=g(2)*g(4); 
m[1,4]=-g(2)*g(3); 
m[1,5]=0; 
m[1,6]=0; 
m[1,7]=0; 
m[1,8]=0; 
m[2,1]=-g(1)*g(2); 
m[2,2]=g(1)^2; 
m[2,3]=-g(2)^2; 
m[2,4]=g(1)*g(2); 
m[2,5]=0; 
m[2,6]=0; 
m[2,7]=0; 
m[2,8]=0; 
m[3,1]=g(3)*g(4); 
m[3,2]=-g(3)^2; 
m[3,3]=g(4)^2; 
m[3,4]=-g(3)*g(4); 
m[3,5]=0; 
m[3,6]=0; 
m[3,7]=0; 
m[3,8]=0; 
m[4,1]=-g(2)*g(3); 
m[4,2]=g(1)*g(3); 
m[4,3]=-g(2)*g(4); 
m[4,4]=g(1)*g(4); 
m[4,5]=0; 
m[4,6]=0; 
m[4,7]=0; 
m[4,8]=0; 
m[5,1]=0; 
m[5,2]=0; 
m[5,3]=0; 
m[5,4]=0; 
m[5,5]=g(1)*g(4); 
m[5,6]=-g(1)*g(3); 
m[5,7]=g(2)*g(4); 
m[5,8]=-g(2)*g(3); 
m[6,1]=0; 
m[6,2]=0; 
m[6,3]=0; 
m[6,4]=0; 
m[6,5]=-g(1)*g(2); 
m[6,6]=g(1)^2; 
m[6,7]=-g(2)^2; 
m[6,8]=g(1)*g(2); 
m[7,1]=0; 
m[7,2]=0; 
m[7,3]=0; 
m[7,4]=0; 
m[7,5]=g(3)*g(4); 
m[7,6]=-g(3)^2; 
m[7,7]=g(4)^2; 
m[7,8]=-g(3)*g(4); 
m[8,1]=0; 
m[8,2]=0; 
m[8,3]=0; 
m[8,4]=0; 
m[8,5]=-g(2)*g(3); 
m[8,6]=g(1)*g(3); 
m[8,7]=-g(2)*g(4); 
m[8,8]=g(1)*g(4); 
*/






