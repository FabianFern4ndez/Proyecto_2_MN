#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#define ERROR 0.000001

void presentacion()
{
   printf("buenos dias\n");	
}

bool comprobadorError(int cant_errores, float errores[]) 
{
    for (int i = 0; i < cant_errores; i++) {
        if (errores[i] > ERROR) {
            return false;
        }
    }
    return true;
}

void jacobi() 
{
    int i, j;
    int filas;
    printf("\nCantidad de variables del sistema: ");
    scanf("%d", &filas);
    int columnas = filas + 1; // Columna adicional para el término independiente

    float matriz[filas][columnas]; 

    // Llenado de la matriz
    for (i = 0; i < filas; i++) {
        printf("Ecuacion %d: \n", i + 1);
        for (j = 0; j < columnas; j++) {
            if (j < filas) {
                printf("Coeficiente de x[%d]: ", j + 1);
                scanf("%f", &matriz[i][j]);
            } else {
                printf("Coeficiente independiente: ");
                scanf("%f", &matriz[i][j]);
            }
        }
    }

    // Mostrar el sistema
    printf("\nSistema de Ecuaciones\n\n");
    for (i = 0; i < filas; i++) {
        for (j = 0; j < columnas; j++) {
            if (j < filas) {
                printf("%fx[%d] ", matriz[i][j], j + 1);
            } else {
                printf("= %f ", matriz[i][j]);
            }
        }
        printf("\n");
    }
    printf("\n");

    // Método de Jacobi
    float errores[filas];
    float x[filas];
    float prev_x[filas];
    for (i = 0; i < filas; i++) {
        errores[i] = 1; // Inicializa errores grandes
        x[i] = 0; // Valores iniciales
    }

    int contador = 0;
    printf("Iteración | ");
    for (i = 0; i < filas; i++) {
        printf("x[%d]       | ", i + 1);
    }
    printf("Errores\n");
    printf("----------------------------------------------------------\n");

    while (!comprobadorError(filas, errores)) {
        // Impresión como en Excel
        printf("%-10d | ", contador);
        for (i = 0; i < filas; i++) {
            printf("%-10f | ", x[i]);
        }

        // Actualización de los valores
        for (i = 0; i < filas; i++) {
            prev_x[i] = x[i]; // Guardar valor anterior
            float suma = 0;

            // Sumar los términos de la ecuación
            for (j = 0; j < filas; j++) {
                if (j != i) {
                    suma += matriz[i][j] * prev_x[j]; // Usa el valor anterior
                }
            }

            // Calcular el nuevo valor
            x[i] = (matriz[i][filas] - suma) / matriz[i][i]; // Término independiente

            // Calcular el error
            errores[i] = fabs(x[i] - prev_x[i]);
            prev_x[i] = x[i]; //profe si está leyendo hasta acá está linea no salía en google y casi
                           //pierdo mi estabilidad mental :'(
        }

        // Mostrar errores
        for (i = 0; i < filas; i++) {
            printf("%-10f | ", errores[i]);
        }
        printf("\n");

        contador++; // Incrementar el contador de iteraciones
    }

    // Mostrar resultados finales
    printf("\nResultados finales:\n");
    for (i = 0; i < filas; i++) {
        printf("x[%d] = %f\n", i + 1, x[i]);
    }
}

void doolittle() {
    int i, j, k;
    int tam_matriz;

    printf("\nCantidad de variables del sistema: ");
    scanf("%d", &tam_matriz);

    // Matriz aumentada: coeficientes + términos independientes
    float matriz[tam_matriz][tam_matriz + 1];

    // Llenado de la matriz
    for (i = 0; i < tam_matriz; i++) {
        printf("Ecuacion %d: \n", i + 1);
        for (j = 0; j < tam_matriz + 1; j++) {
            if (j < tam_matriz) {
                printf("Coeficiente de x[%d]: ", j + 1);
                scanf("%f", &matriz[i][j]);
            } else {
                printf("Coeficiente independiente: ");
                scanf("%f", &matriz[i][j]);
            }
        }
    }

    // Mostrar el sistema
    printf("\nSistema de Ecuaciones\n\n");
    for (i = 0; i < tam_matriz; i++) {
        for (j = 0; j < tam_matriz + 1; j++) {
            if (j < tam_matriz) {
                printf("%fx[%d] ", matriz[i][j], j + 1);
            } else {
                printf("= %f ", matriz[i][j]);
            }
        }
        printf("\n");
    }
    printf("\n");

    // Inicializar L y U
    float L[tam_matriz][tam_matriz];
    float U[tam_matriz][tam_matriz];

    for (i = 0; i < tam_matriz; i++) {
        for (j = 0; j < tam_matriz; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    // Descomposición de Doolittle
    for (i = 0; i < tam_matriz; i++) {
        for (j = i; j < tam_matriz; j++) {
            float sumU = 0;
            for (k = 0; k < i; k++) {
                sumU += L[i][k] * U[k][j];
            }
            U[i][j] = matriz[i][j] - sumU;
        }
        for (j = i + 1; j < tam_matriz; j++) {
            float sumL = 0;
            for (k = 0; k < i; k++) {
                sumL += L[j][k] * U[k][i];
            }
            L[j][i] = (matriz[j][i] - sumL) / U[i][i];
        }
    }

    // Mostrar L y U
    printf("\nMatriz L:\n");
    for (i = 0; i < tam_matriz; i++) {
        for (j = 0; j < tam_matriz; j++) {
            printf("%f ", L[i][j]);
        }
        printf("\n");
    }

    printf("\nMatriz U:\n");
    for (i = 0; i < tam_matriz; i++) {
        for (j = 0; j < tam_matriz; j++) {
            printf("%f ", U[i][j]);
        }
        printf("\n");
    }

    // Resolviendo Ly = b
    float y[tam_matriz];
    for (i = 0; i < tam_matriz; i++) {
        float sum = 0;
        for (j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = matriz[i][tam_matriz] - sum;
    }

    // Resolviendo Ux = y
    float x[tam_matriz];
    for (i = tam_matriz - 1; i >= 0; i--) {
        float sum = 0;
        for (j = i + 1; j < tam_matriz; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    // Mostrar resultados
    printf("\nResultados:\n");
    for (i = 0; i < tam_matriz; i++) {
        printf("x[%d] = %f\n", i + 1, x[i]);
    }
}


void intervaloMedio()
{
	printf("Trabajando en el metodo de Intervalo Medio\n");
}
void regulaFalsi()
{
	printf("Trabajando en el metodo de Regula Falsi\n");
}
float evaluarPolinomio(float polinomio[], int grado, float x)
{
    float resultado = 0.0;
    for (int i = grado; i >= 0; i--) {
        resultado += polinomio[i] * pow(x, i); //esta funcion esta dando problemas pow
    }
    return resultado;
}
float evaluarDerivada(float derivada[], int grado, float x) {
    float resultado = 0.0;
    for (int i = grado - 1; i >= 0; i--) {
        resultado += derivada[i] * pow(x, i);
    }
    return resultado;
}
void newtonRaphson()
{
    int grado;
    int i;
    int x;
	printf("Ingrese el grado del polinomio:");
    scanf("%d",&grado);
    float polinomio[grado+1]; //por el termino lineal
    float derivada[grado];
     //leer el polinomio
    for(i=grado;i>=0;i--)
    {
        printf("Coeficiente de x^%d : ",i);
        scanf("%f" , &polinomio[i]); //bueno se va almacenar el polinomio al reves
    }
    //hallar la derivada por regla de la potencia
    for(i=grado;i>0;i--)
    {
        derivada[i-1] = polinomio[i] * i;
    }
    //imprimir el polinomio
    printf("\nPolinomio: \n");
    for(i=grado;i>=0;i--)
    {
        if(i!=0)
        {
            if(polinomio[i]>=0)
            {
                printf("+%fx^%d  ", polinomio[i], i); 
            }
            else
            {
               printf("%fx^%d  ", polinomio[i], i);  
           }    
        }
        else
        {
           if(polinomio[i]>=0)
            {
                printf("+%f  ", polinomio[i]); 
            }
            else
            {
               printf("%f  ", polinomio[i]);  
           }  
        }
        
    }
    //imprimir derivada
    printf("\nDerivada: \n");
    for(i=grado-1;i>=0;i--)
    {
       if(i!=0)
        {
            if(derivada[i]>=0)
            {
                printf("+%fx^%d  ", derivada[i], i); 
            }
            else
            {
               printf("%fx^%d  ", derivada[i], i);  
           }    
        }
        else
        {
           if(derivada[i]>=0)
            {
                printf("+%f ", derivada[i]); 
            }
            else
            {
               printf("%f ", derivada[i]);  
           }  
        }
    }

    //Metodo de Newton raphson
    float xi;
    printf("\nIngrese el valor inicial xi: ");
    scanf("%f", &xi);

    printf("\nIteración\tXi\t\tf(Xi)\t\tf'(Xi)\t\tError\n");
    int iteracion = 0;
    float error = 1.0;

    do {
        float fxi = evaluarPolinomio(polinomio, grado, xi);
        float fprimaXi = evaluarDerivada(derivada, grado, xi);
        
        if (fprimaXi == 0) {
            printf("La derivada es cero. No se puede continuar.\n");
            return;
        }

        float xiNuevo = xi - fxi / fprimaXi;
        error = fabs(xiNuevo - xi);

        printf("%d\t\t%f\t%f\t%f\t%f\n", iteracion + 1, xi, fxi, fprimaXi, error);

        xi = xiNuevo;
        iteracion++;
        
    } while (error > ERROR && iteracion < 100); // Limitar el número de iteraciones

    printf("Raíz aproximada: %f\n\n", xi);

}
void secante()
{
	printf("Trabajando en el metodo de secante\n");
}
int main()
{
	int decision = 1;
	presentacion();
	while (decision != 0)
	{
		printf("Seleccione un metodo \n");
		printf("1. Metodo de Jacobi\n");
		printf("2. Metodo de Doolittle\n");
        printf("Polinomios\n");
		printf("3. Metodo de Intervalo Medio\n");
		printf("4. Metodo de Regula Falsi\n");
		printf("5. Metodo de Newton Raphson\n");
		printf("6. Metodo de la Secante\n");
		printf("0. Salir del programa\n");
		printf("Eleccion: ");
		scanf("%d", &decision);

		switch (decision)
		{
		case 1:
			jacobi();
			break;
		case 2:
			doolittle();
			break;
		case 3:
			intervaloMedio();
			break;
		case 4:
			regulaFalsi();
			break;
		case 5:
			newtonRaphson();
			break;
		case 6:
			secante();
			break;
		}
	}
}