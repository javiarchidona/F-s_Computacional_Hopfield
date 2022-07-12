//Programa para desordenar un patron aleatoriamente

#include<iostream>
#include<cmath>
#include<fstream>
#include <time.h>
#include <random>  //Libreria que genera numeros aleatorios reales y enteros


#define N 100
#define p 10000        //Porcentaje deformacion de patron (como N=100, p=100%)

using namespace std;

int main(){

    int i, j, k, s[N][N];
    

    ifstream de_i;
    ofstream de_f;


    de_i.open("prueba.txt");
    de_f.open("deformado.txt");



    //Primero leemos el patron original
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            de_i>>s[i][j];
        }
    }

    

    //Generamos posicion aleatoria i, j y la deformamos
    mt19937 semilla(time(NULL));
    uniform_int_distribution<int> random_entero(0,N-1);


    for(k=0;k<p;k++){
        i=random_entero(semilla);
        j=random_entero(semilla);

        s[i][j]=1-s[i][j];
    }



    //Volcamos el patron deformado en su correspondiente archivo
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            de_f<<s[i][j]<<"\t";
        }
        de_f<<endl;
    }
    

    de_i.close();
    de_f.close();


    return 0;
}
