#include<iostream>
#include<cmath>
#include<algorithm>
#include<complex>
#include<cfloat>
#include<iomanip>
#include<vector>
using namespace std;

//

double current =10.0;
double u0=4*(3.14)*pow(10,-7);
double  h=0.005;
double poles=2;
double slots = 12;
double k_w=0.955;
double conductor_per_slot = 100;
double len_of_primary = 0.15;
double pole_pitch = len_of_primary/poles;
double k=3.14*pole_pitch;
double f=50;
double v=3;
double omega=2*3.14*f;
double q=0.2;
double sigma=2.72*pow(10,7);
double s1=0.001;
double s2=1;
double s3=0.001;
double s4=1000;
double j_max=(sqrt(2)*current*slots*conductor_per_slot*k_w)/len_of_primary;
complex<double> Gamma(2+2*q*q,omega*u0*sigma*q*q*h*h);
double Alpha=(q*q)*(1+0.5*u0*sigma*v*h);
double Beta=q*q*(1-0.5*u0*sigma*v*h);
double t=0;
int point=10;

//


void display(vector<vector<complex< double>>>v, int row,int col){
    for( int i{0};i<row;i++){
        cout<<"J = "<<46-i<<" : "<<endl;
        for( int j{0};j<col;j++){
        cout<<v[i][j];}
        cout<<endl;
        cout<<endl;
        cout<<endl;
    }
    cout<<endl;
}

 complex<double> laplace(vector<vector<complex<double>>>v, int row,int col){
    if(row==0){
        if(col==0)
             return ((v[row][col+1]+polar(q*q,0.0)*(v[row+1][col]))/(2*(1+q*q)));
        else if(col==37)
             return ((v[row][col-1]+polar(q*q,0.0)*(v[row+1][col]))/(2*(1+q*q)));
        else
             return ((v[row][col-1]+v[row][col+1]+polar(q*q,0.0)*(v[row+1][col]))/(2*(1+q*q)));
    }
    else if(col=0)
         return ((v[row][col+1]+polar(q*q,0.0)*(v[row-1][col]+v[row+1][col]))/(2*(1+q*q)));
    else if(col=37)
         return ((v[row][col-1]+polar(q*q,0.0)*(v[row-1][col]+v[row+1][col]))/(2*(1+q*q)));
    else
        return ((v[row][col-1]+v[row][col+1]+polar(q*q,0.0)*(v[row-1][col]+v[row+1][col]))/(2*(1+q*q)));
 }

 complex<double> interference_iron_al(vector<vector<complex< double>>>v, int row,int col){
    if(col==0)
        return (((Beta+s1*q*q)*v[row][col+1]+polar(2.0,0.0)*v[row-1][col]+polar(2.0,0.0)*v[row+1][col])/(Gamma+2*s1*(1+q*q)));
    else if(col==37)
        return (((Alpha+s1*q*q)*(v[row][col-1])+polar(2.0,0.0)*v[row-1][col]+polar(2.0,0.0)*v[row+1][col])/(Gamma+2*s1*(1+q*q)));
    else
        return (((Alpha+s1*q*q)*(v[row][col-1])+(Beta+s1*q*q)*v[row][col+1]+polar(2.0,0.0)*v[row-1][col]+polar(2.0,0.0)*v[row+1][col])/(Gamma+2*s1*(1+q*q)));
}

 complex<double> al_poisson(vector<vector<complex< double>>>v, int row,int col){
    if(col==0)
        ((polar(Beta,0.0)*v[row][col+1]+v[row-1][col]+v[row+1][col])/(Gamma+2*s1*(1+q*q)));
    else if(col==37)
        return ((polar(Alpha,0.0)*v[row][col-1]+v[row-1][col]+v[row+1][col])/(Gamma+2*s1*(1+q*q)));
    else
        return ((polar(Alpha,0.0)*v[row][col-1]+polar(Beta,0.0)*v[row][col+1]+v[row-1][col]+v[row+1][col])/(Gamma+2*s1*(1+q*q)));
}

 complex<double> al_air_interference(vector<vector<complex< double>>>v,int row,int col){
    if(col==0)
        return (((Beta*s2+q*q)*v[row][col+1]+polar(2.0,0.0)*v[row-1][col]+polar(2.0,0.0)*v[row+1][col])/(Gamma*s2+2*(1+q*q)));
    else if(col==37)
        return (((Alpha*s2+q*q)*v[row][col-1]+polar(2.0,0.0)*v[row-1][col]+polar(2.0,0.0)*v[row+1][col])/(Gamma*s2+2*(1+q*q)));
    else
        return (((Alpha*s2+q*q)*v[row][col-1]+(Beta*s2+q*q)*v[row][col+1]+polar(2.0,0.0)*v[row-1][col]+polar(2.0,0.0)*v[row+1][col])/(Gamma*s2+2*(1+q*q)));
}

 complex<double> current_sheet(vector<vector<complex< double>>>v,int row,int col){
        return ((s3*v[row-1][col]+v[row+1][col]-u0*q*h*j_max*cos(omega*t-k*h))/(s3+1));
}

 complex<double> corner(vector<vector<complex< double>>>v,int row,int col){
    if(col==0)
        return ((((v[row][col+1])+(s3+1)*(v[row-1][col]+v[row+1][col]))/(4*(s3+1)*(1+q*q)))+((s3*v[row-1][col]+v[row+1][col]-u0*q*h*cos(omega*t-k*h))/(2*s3+2)));
    else if(col==37)
        return (((2*q*q*(s3*v[row][col-1])+(s3+1)*(v[row-1][col]+v[row+1][col]))/(4*(s3+1)*(1+q*q)))+(s3*v[row-1][col]+v[row+1][col]-u0*q*h*cos(omega*t-k*h))/(2*s3+2));
    else
        return (((2*q*q*(s3*v[row][col-1]+v[row][col+1])+(s3+1)*(v[row-1][col]+v[row+1][col]))/(4*(s3+1)*(1+q*q)))+(s3*v[row-1][col]+v[row+1][col]-u0*q*h*cos(omega*t-k*h))/(2*s3+2));
}

complex<double> LHS_air_and_core(vector<vector<complex< double>>>v,int row,int col){
        return (((2*q*q)*(s4*v[row][col-1]+v[row][col+1])+(s4+1)*(v[row+1][col]))/((s4+1)*2*(1+q*q)));
}

 complex<double> LHS_core_and_air(vector<vector<complex< double>>>v,int row,int col){
        return (((2*q*q)*(s3*v[row][col-1]+v[row][col+1])+(s3+1)*(v[row+1][col]))/((s3+1)*2*(1+q*q)));
}

complex<double> Bx(vector<vector<complex<double>>>v, int row,int col){
    
    if(col=0)
         return ((v[row][col+1])/(2*q*h));
    else if(col=37)
         return ((-v[row][col-1])/(2*q*h));
    else
        return ((v[row][col+1]-v[row][col-1])/(2*q*h));
 }

complex<double> By(vector<vector<complex<double>>>v, int row,int col){
        return ((v[row-1][col]-v[row+1][col])/(2*h));
 }



//

 int main()
{
	  vector<vector<complex< double>>> v( 46 , vector<complex< double>> (38, 0));
       vector<vector<complex< double>>> util( 46 , vector<complex< double>> (38, 0));
        vector<complex< double>>residue(100,0);
      vector<complex< double>>vec_Bx (38, 0);
       vector<complex< double>>vec_By (38, 0);



     char z;
while(true){
  cout<<"press a/A for t=0 or b/B for t=T/4 : ";
  cin>>z;
  if(z=='a'|| z=='A'){
  t=0;
  break;}
  else if(z=='b'|| z=='B'){
  t=1/(4*f);
  break;}
  else{
      cout<<"Try again!"<<endl;}
};


      while(true){
       int temp;
      cout<<"Enter point of reference for which  we store calculated residue (between 1 and 37) : ";
      cin>>temp;
      if(temp>=1 && temp<=37){
          point=temp;
          break;
      }
      else
          cout<<"try again !"<<endl;
      }
      
      int k=1;
      while(k<100){


       int row,col;

      //for J=1 :not req.

      //for J=2,3
      for(row=44;row>=43;row--){
          for(col=0;col<38;col++){
              v[row][col]=laplace(v,row,col);
          }
      }
      //for J=4
      row=42;
     for(col=0;col<38;col++)
              v[row][col]=interference_iron_al(v,row,col);

      //J=5,6
      for(row=41;row>=40;row--){
          for(col=0;col<38;col++)
              v[row][col]=al_poisson(v,row,col);
      }

      //J=7
      row=39;
       for(col=0;col<38;col++)
              v[row][col]=al_air_interference(v,row,col);

      //for J=8 to J=15
      for(row=38;row>=31;row--){
          for(col=0;col<38;col++){
              v[row][col]=laplace(v,row,col);
          }
      }

      //for J=16
      row=30;
          for(col=0;col<=9;col++){
              v[row][col]=laplace(v,row,col);
          }
          v[row][10]=corner(v,row,col);
          for(col=11;col<=25;col++)
              v[row][col]=current_sheet(v,row,col);
            v[row][26]=corner(v,row,col);
            for(col=27;col<38;col++)
                v[row][col]=laplace(v,row,col);
            for(col=0;col<38;col++){
              vec_Bx[col]=Bx(v,row,col);
              vec_By[col]=By(v,row,col);
          }



          //for J=17 to 46
      for(row=29;row>=0;row--){
          for(col=0;col<=9;col++){
              v[row][col]=laplace(v,row,col);
          }
          v[row][10]=LHS_air_and_core(v,row,col);
          for(col=11;col<=25;col++)
              v[row][col]=laplace(v,row,col);
            v[row][26]=LHS_core_and_air(v,row,col);
            for(col=27;col<38;col++)
                v[row][col]=laplace(v,row,col);
      }


      
      
      //Updated value = calculate-alfa*residue
      for(row=44;row>=0;row--){
          for(col=0;col<38;col++){
              v[row][col]=v[row][col]-polar(1.4,-0.03)*(v[row][col]-util[row][col]);
          }
      }
      //residue=calculated-old
     row=30;col=point;
      residue[k]=v[row][col]-util[row][col];
      
      //store for next
      for(row=44;row>=0;row--){
          for(col=0;col<38;col++){
             util[row][col]=v[row][col];
          }
      }
      
      k++;
      }
      display(v,46,38);
      

      cout<<"Abs(Az) : "<<endl;
      int row,col;
      row=30;
          for(col=0;col<38;col++){
              cout<<abs(v[row][col])<<"  ";
          }
      
      cout<<endl;
      cout<<endl;
      cout<<endl;
      
      cout<<"Bx :"<<endl;
      for(auto i: vec_Bx)
          cout<<i;
    cout<<endl;
    cout<<endl;
    cout<<"Abs(Bx) :"<<endl;
    for(auto i: vec_Bx)
          cout<<abs(i)<<" ";
          
    cout<<endl;
    cout<<endl;
    cout<<"By :"<<endl;
      for(auto i: vec_By)
          cout<<i;
    cout<<endl;
    cout<<endl;
    cout<<"Abs(By) :"<<endl;
    for(auto i: vec_By)
          cout<<abs(i)<<" ";
      cout<<endl;
      cout<<endl;
      cout<<"Residue at J=16 and point "<<point<<" :"<<endl;
     for(auto i:residue){
         cout<<i;
     }
         
         cout<<endl;
         cout<<endl;
      
      
      return 0;
}