LIB("rinvar.lib"); 
LIB("invar.lib"); 
LIB("matrix.lib");

int N = 2; 


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
		A2[i,j] = a(N*N+(i-1)*N+j); 
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







