#include <stdio.h>
#include <math.h>
#include <windows.h>

typedef struct
{
    double a;
    double b;
    double fa;
    double fb;
    double p;
    double fp;
    double er;
}biseccion;

typedef struct
{
    double x;
    double fx;
    double dx;
    double er;
}newton;

double fx1(double x)
{
    return(pow(x,2)*cos(x)-2*x);
}
double dx1(double x)
{
    return(2*x*cos(x)-pow(x,2)*sin(x)-2);
}
double fx2(double x)
{
    double e = 2.71828;
    return((6-(2/pow(x,2)))*((pow(e,2+x))/4)+1);
}
double dx2(double x)
{
    double e = 2.71828;
    return(((pow(e,x+2))*((3*pow(x,3))-x+2))/(2*pow(x,3)));
}
double fx3(double x)
{
    return(pow(x,3) - 3*sin(pow(x,2))+1);
}
double dx3(double x)
{
    return(3*pow(x,2)-6*x*cos(pow(x,2)));
}
double fx4(double x)
{
    return((pow(x, 3)) + (6*pow(x,2)) + (9.4*x) + 2.5);
}
double dx4(double x)
{
    return((3 * pow(x, 2)) + (12*x) + 9.4);
}

void impN(newton funcion, int i)
{
    printf("Iteracion %d: x = %lf, fx = %lf, dx = %lf, Error relativo = %lf\n", i+1, funcion.x, funcion.fx, funcion.dx, funcion.er);
}

void impB(biseccion funcion, int i)
{
    printf("Iteracion %d: a = %lf, b = %lf, f(a) = %lf, f(b) = %lf, p = %lf, f(p) = %lf, Error relav = %lf\n", i+1, funcion.a, funcion.b, funcion.fa, funcion.fb, funcion.p, funcion.fp, funcion.er);
}

double absoluto(double v)
{
    if(v < 0)
    {
        v *= -1;
    }
    return(v);
}

double errorR(double x0, double x1)
{
    double p1 = (absoluto(x1 - x0));
    double p2 = absoluto(x1);
    return(p1/p2);
}

int metodoBiseccion(int f)
{
    biseccion funcion;
    int iteraciones = 0;
    double a = 0, b = 0;
    double tolerancia = 0;
    double aux = 0;
    double ra = 0;
    printf("Ingrese numero de iteraciones: ");
    scanf("%d", &iteraciones);
    do
    {
        printf("Ingrese su intervalo, inicio: ");
        scanf("%lf", &funcion.a);
        printf("a fin: ");
        scanf("%lf", &funcion.b);
        switch(f)
        {
            case 1:
                a = fx1(funcion.a);
                b = fx1(funcion.b);
            break;

            case 2:
                a = fx2(funcion.a);
                b = fx2(funcion.b);
            break;

            case 3:
                a = fx3(funcion.a);
                b = fx3(funcion.b);
            break;

            case 4:
                a = fx4(funcion.a);
                b = fx4(funcion.b);
            break;

            default:
            break;
        }
        if(a*b> 0)
        {
            printf("La raiz no esta en este intervalo, ingrese un nuevo intervalo\n");
        }
    }while(a*b> 0);
    printf("Ingrese tolerancia: ");
    scanf("%lf", &tolerancia);
    for(int i = 0; i < iteraciones; i++)
    {
        funcion.p = (funcion.a + funcion.b)/2;
        switch(f)
        {
            case 1:
                funcion.fa = fx1(funcion.a);
                funcion.fb = fx1(funcion.b);
                funcion.fp = fx1(funcion.p);
            break;

            case 2:
                funcion.fa = fx2(funcion.a);
                funcion.fb = fx2(funcion.b);
                funcion.fp = fx2(funcion.p);
            break;

            case 3:
                funcion.fa = fx3(funcion.a);
                funcion.fb = fx3(funcion.b);
                funcion.fp = fx3(funcion.p);
            break;

            case 4:
                funcion.fa = fx4(funcion.a);
                funcion.fb = fx4(funcion.b);
                funcion.fp = fx4(funcion.p);
            break;

            default:
            break;
        }
        funcion.er = errorR(aux, funcion.p);
        aux = funcion.p;
        impB(funcion, i);
        if(i > 0)
        {
            if(funcion.er < tolerancia)
            {
                printf("Se revaso la tolerancia\n");
                printf("La aproximacion de raiz encontrada es: %lf en la iteracion %d", aux, i+1);
                break;
            }
            else if(i == iteraciones-1)
            {
                printf("No se revaso la tolerancia\n");
                printf("La aproximacion de raiz encontrada es: %lf, en la iteracion %d\n", aux, i+1);
                break;
            }
        }
        if(funcion.fa*funcion.fp < 0)
        {
            funcion.b = funcion.p;
        }
        else
        {
            funcion.a = funcion.p;
        }
    }
}

int metodoNewton(int f)
{
    int iteraciones = 0;
    double aux = 0;
    double tolerancia = 0;
    newton funcion;
    printf("Ingrese el numero de iteraciones que desea hacer\n");
    scanf("%d", &iteraciones);
    printf("Ingrese el valor inicial\n");
    scanf("%lf", &funcion.x);
    printf("Ingrese la tolerancia\n");
    scanf("%lf", &tolerancia);
    funcion.er = 0;
    for(int i = 0; i < iteraciones; i++)
    {
        switch(f)
        {
            case 1:
                funcion.fx = fx1(funcion.x);
                funcion.dx = dx1(funcion.x);
            break;

            case 2:
                funcion.fx = fx2(funcion.x);
                funcion.dx = dx2(funcion.x);
            break;

            case 3:
                funcion.fx = fx3(funcion.x);
                funcion.dx = dx3(funcion.x);
            break;

            case 4:
                funcion.fx = fx4(funcion.x);
                funcion.dx = dx4(funcion.x);
            break;

            default:
                break;
        }
        impN(funcion, i);
        aux = funcion.x;
        funcion.x = funcion.x-(funcion.fx/funcion.dx);
        funcion.er = errorR(aux, funcion.x);
        if(i > 0)
        {
            if(funcion.er < tolerancia)
            {
                printf("Se revaso la tolerancia\n");
                printf("La aproximacion de raiz encontrada es: %lf, en la iteracion %d\n", aux, i+1);
                break;
            }
            else if(i == iteraciones-1)
            {
                printf("No se revaso la tolerancia\n");
                printf("La aproximacion de raiz encontrada es: %lf, en la iteracion %d\n", aux, i+1);
                break;
            }
        }
    }
}

int d;

void imprimir(float matriz[d][d])
{
    printf("\n");
    for(int x = 0; x < d; x++)
    {
        printf("| ");
        for(int y = 0; y < d; y++)
        {
            printf("%f ", matriz[x][y]);
        }
        printf(" |\n");
    }
}

void imprimirMV(float matriz[d][d], float vector[d])
{
    printf("\n");
    for(int x = 0; x < d; x++)
    {
        printf("| ");
        for(int y = 0; y < d; y++)
        {
            printf("%f ", matriz[x][y]);
        }
        printf("|%f|\n", vector[x]);
    }
}

void triangularSup(float matriz[d][d])
{
    int x = 0, y = 0, z = 0;
    float aux = 0, aux2 = 0;
    for(z = 0; z < d; z++)
    {
        for(x = 0; x < d; x++)
        {
            aux = matriz[x][z];
            for(y = 0; y < d; y++)
            {
                if(x > z)
                {
                    aux2 = aux*matriz[z][y]/matriz[z][z];
                    matriz[x][y] = matriz[x][y] - aux2;
                }
            }
        }
    }
    imprimir(matriz);
}
int determinante(float m[d][d]){

    triangularSup(m);
    float det = 1;
    for(int x = 0; x < d; x++)
    {
        for(int y = 0; y < d; y++)
        {
            if(x == y)
            {
                det *= m[x][y];
            }
        }
    }
    printf("Det: %.2f\n", det);
    return det;
}

int dimensionMatriz()
{
    int v = 0;
    do
    {
        printf("Ingresa la dimension de la matriz\n");
        scanf("%d", &v);

        if(v <= 0)
        {
            printf("Valor incorrecto, favor de ingresar nuevamente\n");
        }
    }while(v <= 0);
    return v;
}

void leeM(float matriz[d][d])
{
    char opc;
    int c1 = 0, c2 = 0;
    for(int x = 0; x < d; x++)
    {
        for(int y = 0; y < d; y++)
        {
            printf("Ingresa [%d][%d]: ", x+1, y+1);
            scanf("%f", &matriz[x][y]);
        }
    }
    imprimir(matriz);
    fflush(stdin);
    printf("\nDesas corregir algun dato de su matriz? s -> SI, cualquier tecla -> NO\n");
    scanf("%c", &opc);
    while(opc == 's')
    {
        printf("Ingresa las coordenadas del dato que quieres corregir\n");
        do
        {
            printf("Numero de renglon: ");
            scanf("%d", &c1);
            if(c1 < 1 || c1 > d)
            {
                printf("Numero invalido, ingresa nuevamente\n");
            }
        }while(c1 < 1 || c1 > d);
        do
        {
            printf("Numero de columna: ");
            scanf("%d", &c2);
            if(c2 < 1 || c2 > d)
            {
                printf("Numero invalido, ingresa nuevamente\n");
            }
        }while(c2 < 1 || c2 > d);

        printf("Por que valor deseas cambiar el valor '%.2f': ", matriz[c1-1][c2-1]);
        scanf("%f", &matriz[c1-1][c2-1]);

        imprimir(matriz);

        fflush(stdin);
        printf("\nDesas corregir algun dato? s -> SI, cualquier tecla -> NO\n");
        scanf("%c", &opc);
    }
}

void imprimirV(float vector[d]){

    printf("\n");
    for(int x = 0; x < d; x++){
        printf("| %f |\n", vector[x]);
    }
    printf("\n");
}

void leeV(float vector[d]){

    char opc;
    int c = 0;
    for(int x = 0; x < d; x++){
        printf("Ingrese el valor de la posicion [%d] de su vector: ", x+1);
        scanf("%f", &vector[x]);
    }
    imprimirV(vector);
    fflush(stdin);
    printf("Deseas corregir algun dato de su vector? s -> SI, cualquier tecla -> NO\n");
    scanf("%c", &opc);
    while(opc == 's'){
        printf("Ingresa las coordena del valor que quieres corregir: ");
        do
        {
            scanf("%d", &c);

            if(c < 1 || c > d)
            {
                printf("Valor incorrecto, ingrese nuevamente\n");
            }
        }while(c < 1 || c > d);
        printf("Por que valor deseas cambiar el valor '%.2f': ", vector[c-1]);
        scanf("%f", &vector[c-1]);

        printf("Cambio exitoso\n");

        imprimirV(vector);

        fflush(stdin);
        printf("Deseas corregir algun dato de su vector? s -> SI, cualquier tecla -> NO\n");
        scanf("%c", &opc);
    }

}

void verificaEDD(float matriz[d][d])
{
    int aux = 0;
    int suma = 0;

    for(int x = 0; x < d; x++)
    {
        for(int y = 0; y < d; y++)
        {
            if(x != y)
            {
                suma += matriz[x][y];
            }
        }
        if(suma > matriz[x][x])
        {
            printf("La convergencia no se garantiza por no tratarse de un sistema EDD\n");
            break;
        }else{
            suma = 0;
        }
    }
    printf("EDD\n");
}

void programa1()
{
    int met = 0;
    int fun = 0;
    printf("Hola, selecciona la funcion que deseas utilizar\n");

    do
    {
        printf("1.- f(x) = x^2cos(x)-2x\n");
        printf("2.- f(x) = (6-(2/x^2))*(e^(2+x)/4)+1\n");
        printf("3.- f(x) = x^3-3*sen(x^2)+1\n");
        printf("4.- f(x) = x^3+6x^2+9.4x+2.5\n");
        printf("5.- Salir\n");
        scanf("%d", &fun);
        if(fun > 0 && fun < 5)
        {
            do
            {
                printf("\n1.- Metodo de Biseccion\n");
                printf("2.- Metodo de Newton\n");
                printf("3.- Limpiar la pantalla\n");
                printf("4.- Regresar a funciones\n");
                scanf("%d", &met);
                switch(met)
                {
                    case 1:
                        metodoBiseccion(fun);
                    break;

                    case 2:
                        metodoNewton(fun);
                    break;

                    case 3:
                        system("cls");
                    break;

                    case 4:
                        printf("Escoga una nueva funcion\n");
                    break;
                }
            }while(met != 4);
        }
    }while(fun != 5);
}

void despeje(float matriz[d][d], float vector[d])
{
    int aux = 0;
    for(int x = 0; x < d; x++)
    {
        aux = matriz[x][x];
        for(int y = 0; y < d; y++)
        {
            if(x == y)
            {
                matriz[x][y] = 0;
            }
            else{
                matriz[x][y] /= -aux;
            }
        }
        vector[x] /= aux;
    }
    imprimirMV(matriz, vector);
}

float normaEsp(float vector[d])
{
    float mayor = 0;
    for(int x = 0; x < d; x++)
    {
        if(vector[x] < 0)
        {
            vector[x] *= -1;
        }
    }
    mayor = vector[0];
    for(int x = 0; x < d; x++){
        if(vector[x] > mayor)
        {
            mayor = vector[x];
        }
    }
    printf("La norma espectral es: %f\n", mayor);
    return mayor;
}

float error(float vectorI[d], float vectorN[d])
{
    float vectorE[d];

    for(int x = 0; x < d; x++)
    {
        vectorE[x] = vectorN[x] - vectorI[x];
    }
    printf("Error de la iteracion: \n");

    imprimirV(vectorE);
    return(normaEsp(vectorE));
}

void multiplicaMV(float matriz[d][d], float vector[d], float vectorI[d], int iteraciones, int cont, float tolerancia)
{
    float vectorN[d];
    float normaE;
    float aux = 0;
    for(int x = 0; x < d; x++)
    {
        for(int y = 0; y < d; y++)
        {
            aux += matriz[x][y] * vectorI[y];
        }
        aux += vector[x];
        vectorN[x] = aux;
        aux = 0;
    }
    --iteraciones;
    printf("\nIteracion %d:\n", ++cont);
    imprimirV(vectorN);

    normaE = error(vectorI, vectorN);

    if(iteraciones != 0 && normaE > tolerancia)
    {
        multiplicaMV(matriz, vector, vectorN, iteraciones, cont, tolerancia);
    }else if(iteraciones == 0)
    {
        printf("Se han cumplido todas las iteraciones, el resultado aproximado es:\n");
        imprimirV(vectorN);
    }else{
        printf("Se ha revasado la tolerancia en la iteracion %d, el resultado aproximado es:\n", cont);
        imprimirV(vectorN);
    }
}

void programa3(float matriz[d][d], float vector[d]){ //Jacobi

    float vectori[d];
    int iteraciones = 0;
    float tolerancia = 0.0;
    char opc;
    verificaEDD(matriz);
    despeje(matriz, vector);
    do
    {
        printf("Ingresa tu vector inicial para las iteraciones: \n");
        leeV(vectori);

        printf("Ingresa el numero de iteraciones: ");
        scanf("%d", &iteraciones);

        printf("Ingresa la tolerancia: ");
        scanf("%f", &tolerancia);

        multiplicaMV(matriz, vector, vectori, iteraciones, 0, tolerancia);

        printf("Deseas iterar con otro nuevo vector inicial? pulse 's' para 'SI', pulse cualquier tecla para 'NO': ");
        fflush(stdin);
        scanf("%c", &opc);
    }while(opc == 's');
}

void programa2()
{
    d = dimensionMatriz();
    float matriz[d][d];
    float maux[d][d];
    float vector[d];
    int det = 0;

    leeM(matriz);
    leeV(vector);

    for(int x = 0; x < d; x++)
    {
        for(int y = 0; y < d; y++)
        {
            maux[x][y] = matriz[x][y];
        }
    }

    printf("Diagonal\n");
    det = determinante(maux);

    printf("Matriz principal\n");
    imprimirMV(matriz, vector);


    if(det == 0)
    {
        printf("\nEl determinante es 0, no pueden tiene solucion\n");
    }
    else
    {
        programa3(matriz, vector);
    }
}

int main(){

    int opc = 0;

    printf("Hola, bienvenido a nuestro programa 3!\n");

    do
    {
        printf("Que opcion deseas usar?\n");
        printf("1.- Metodos numericos (biseccion y newton)\n");
        printf("2.- Solucion de sistema de ecuaciones\n");
        printf("3.- Salir\n");
        scanf("%d", &opc);
        switch(opc){
            case 1:
                programa1();
            break;

            case 2:
                programa2();
            break;

            case 3:
                printf("Gracias por usar el programa!\n");
            break;

            default:
                printf("Opcion no valida, favor de ingrese nuevamente\n");
            break;
        }
    }while(opc != 3);

    printf("Elaborado por:\n");
    printf("Kevin Sánchez Sánchez\n");
    printf("Alison Abril Navarro Pérez\n");
    printf("José Carlos Carbajal Mejía")

}
