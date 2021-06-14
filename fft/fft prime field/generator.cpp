#include<cstdio>
#include<iostream>
using namespace std;
const int MAXN=5e6+5;
const int mo=998244353;
const int g=3;
int n,m,limit,len,invn;
int a[MAXN],b[MAXN],p[MAXN];
inline int read()
{
    int x=0,f=1; char ch=getchar(); 
    while(ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
    while(ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
    return x*f;
}
int power(int a,int b)
{
    int res=1;
    while (b)
    {
        if (b&1) res=1ll*res*a%mo;
        a=1ll*a*a%mo;
        b>>=1;
    }
    return res;
}
void NTT(int *A,int type)
{
    for (int i=0; i<limit; i++) 
        if (i<p[i]) swap(A[i],A[p[i]]);
    for (int l=1; l<limit; l<<=1)
    {
        int wn=power(g,(mo-1)/(l<<1));
        if (type==-1) wn=power(wn,mo-2);
        for (int i=0; i<limit; i+=(l<<1))
        {
            int w=1;
            for (int j=0; j<l; j++,w=1ll*w*wn%mo)
            {
                int t=1ll*w*A[i+j+l]%mo;
                A[i+j+l]=(A[i+j]-t+mo)%mo;
                A[i+j]=(A[i+j]+t)%mo;
            }
        }
    }
}
int main()
{
    n=read(); m=read(); 
    for (int i=0; i<=n; i++) a[i]=read();
    for (int i=0; i<=m; i++) b[i]=read();
    
    limit=1;
    while (limit<n+m+1) limit<<=1,len++;
    for (int i=0; i<limit; i++)
        p[i]=(p[i>>1]>>1)|((i&1)<<(len-1));
    
    NTT(a,1); NTT(b,1);
    for (int i=0; i<limit; i++) a[i]=1ll*a[i]*b[i]%mo;
    NTT(a,-1);
    invn=power(limit,mo-2);
    for (int i=0; i<=n+m; i++)
        printf("%lld ",1ll*a[i]*invn%mo);
    printf("\n");
    return 0;
}