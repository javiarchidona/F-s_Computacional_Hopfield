#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <random>  //Libreria que genera numeros aleatorios reales y enteros


#define T 0.1        //Valor temperatura
#define N 100         //Dimension de la matriz 
#define PMC 30        //Pasos Monte Carlo

double omega[N][N][N][N];
double CalcularP(int n, int m, double omega[N][N][N][N], double theta[N][N], double s[N][N]);    //Funcion que calcula el valor de P 



using namespace std;

int main(){

    int i, j, k, l, n, m, z, xi[N][N];
    double random, p, a, sum, chi, sol, num;
    double s[N][N], theta[N][N];


    ofstream datos, e_inicial, e_final, solapamiento,da;
    datos.open("datos.txt");
    e_inicial.open("estado_inicial.txt");
    e_final.open("estado_final.txt");
    solapamiento.open("solapamiento.txt");


    //Patrones a leer
    ifstream pr, de;
    pr.open("prueba.txt");
    de.open("deformado.txt");


    mt19937 semilla(time(NULL));
    uniform_int_distribution<int> random_entero(0,N-1);
    uniform_real_distribution<double> random_real(0.0,1.0);



    //Leemos el patron original y deformado



    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            pr>>xi[i][j];
        }
    }



    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            de>>s[i][j];
        }
    }



    //Se calcula a
    a=0.0;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            a=a+xi[i][j];
        }
    }

    a=a/(N*N);


    //Se calcula omega
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            for(k=0;k<N;k++){
                for(l=0;l<N;l++){
                    if((i==k)&&(j==l)){
                        omega[i][j][k][l]=0.0;
                    }
                    else{
                        omega[i][j][k][l]=(xi[i][j]-a)*(xi[k][l]-a);
                        omega[i][j][k][l]=omega[i][j][k][l]/(N*N);
                    }
                }
            }
        }
    }

    //Se calcula theta
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            sum=0.0;

            for(k=0;k<N;k++){
                for(l=0;l<N;l++){
                    sum=sum+omega[i][j][k][l];
                }
            }
            theta[i][j]=0.5*sum;
            
        }
    }

/*
    //Si queremos partir de un estado inicial aleatorio generamos la matriz s aleatoriamente
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            ran=random_real(semilla);
            if(ran<=0.5)
                s[i][j]=1.0;
            else
                s[i][j]=0.0;
        }
    }

*/


    //Volcamos los datos iniciales en el fichero

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            e_inicial<<s[i][j]<<"\t";
        }
        e_inicial<<endl;
    }






    for(z=0;z<PMC;z++){

        //Imprimimos los datos en archivo .txt

        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                datos<<s[i][j]<<"\t";
            }
            datos<<endl;
        }

        datos<<endl;


        for(i=0;i<(N*N);i++){

            //Elegimos punto (n,m) aleatorio de la red

            n=random_entero(semilla);
            m=random_entero(semilla);

            //Evaluamos p=min[1,exp-(E/T)]

            p=CalcularP(n, m, omega, theta, s);

            //Generamos numero aleatorio uniforme en el intervalo [0,1], y si es <=p cambiamos el signo al elemento

            chi=random_real(semilla);

            if(chi<p){
                if(s[n][m]==1.0){
                    s[n][m]=0.0;
                }
                else    s[n][m]=1.0;
            }
        }


        //Se calcula el solapamiento

        sum=0.0;

        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                sum=sum+(xi[i][j]-a)*(s[i][j]-a);
            }
        }
        sol=sum/(N*N*a*(1.0-a));

        solapamiento<<z<<"\t"<<sol<<endl;

        

    }

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            e_final<<s[i][j]<<"\t";
        }
        e_final<<endl;
    }

 



    datos.close();
    e_inicial.close();
    e_final.close();
    solapamiento.close();
    pr.close();
    de.close();


    return 0;
}


double CalcularP(int n, int m, double omega[N][N][N][N], double theta[N][N], double s[N][N]){

    int i, j, k, l;
    double p, sum, H;

    //Se calcula el hamiltoniano

    sum=0.0;

    for(k=0;k<N;k++){
        for(l=0;l<N;l++){

            if((k!=n)&&(l!=m)){
                sum=sum-(omega[n][m][k][l]*s[k][l]);
            }   

        }
    }


    H=(1-2*s[n][m])*(theta[n][m]+0.5*sum);


    p=min(1.0,exp(-H/T));


    return p;

}
