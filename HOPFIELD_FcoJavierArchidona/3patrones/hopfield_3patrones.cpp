#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <random>  //Libreria que genera numeros aleatorios reales y enteros


#define T 0.0001        //Valor temperatura
#define N 50         //Dimension de la matriz 
#define PMC 20        //Pasos Monte Carlo

double omega[N][N][N][N];
double CalcularP(int n, int m, double omega[N][N][N][N], double theta[N][N], double s[N][N]);    //Funcion que calcula el valor de P 



using namespace std;

int main(){

    int i, j, k, l, n, m, z, xi[N][N][3], mu;
    double random, p, a[3], sum, chi, sol[3], num, ran;
    double s[N][N],  theta[N][N];

    //Archivos a escribir

    ofstream e_inicial, e_final, solapamiento, datos;
    datos.open("datos.txt");
    e_inicial.open("estado_inicial.txt");
    e_final.open("estado_final.txt");
    solapamiento.open("solapamiento.txt");


    //Archivos a leer

    ifstream p1, p2, p3, de;
    p1.open("prueba1.txt");
    p2.open("prueba2.txt");
    p3.open("prueba3.txt");
    de.open("deformado.txt");


    mt19937 semilla(time(NULL));
    uniform_int_distribution<int> random_entero(0,N-1);
    uniform_real_distribution<double> random_real(0.0,1.0);



    //Leemos los tres patrones prueba

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            p1>>xi[i][j][0];
        }
    }


    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            p2>>xi[i][j][1];
        }
    }


    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            p3>>xi[i][j][2];
        }
    }



    //Leemos el patron deformado

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            de>>s[i][j];
        }
    }



    //Se calcula a
    
    for(mu=0;mu<3;mu++){
        a[mu]=0.0;
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                a[mu]=a[mu]+xi[i][j][mu];
            }
        }

        a[mu]=a[mu]/(N*N);
    }



    //Se calcula omega

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            for(k=0;k<N;k++){
                for(l=0;l<N;l++){
                    if((i==k)&&(j==l)){
                        omega[i][j][k][l]=0.0;
                    }
                    else{
                        omega[i][j][k][l]=(xi[i][j][0]-a[0])*(xi[k][l][0]-a[0])+(xi[i][j][1]-a[1])*(xi[k][l][1]-a[1])+(xi[i][j][2]-a[2])*(xi[k][l][2]-a[2]);
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
//Si se quiere empezar desde estado inicial aleatorio se genera la matriz s de forma aleatoria

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




    //Volcamos los datos iniciales en el fichero "estado_inicial.txt"

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

        for(mu=0;mu<3;mu++){

            sum=0.0;

            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    sum=sum+(xi[i][j][mu]-a[mu])*(s[i][j]-a[mu]);
                }
            }

            sol[mu]=sum/(N*N*a[mu]*(1.0-a[mu]));

        }


        solapamiento<<z<<"\t"<<sol[0]<<"\t"<<sol[1]<<"\t"<<sol[2]<<endl;


    }


    //Escribimos el estado final del patron en el archivo "estado_final.txt"

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
    p1.close();
    p2.close();
    p3.close();
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