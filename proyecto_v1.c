#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#define ERROR 0.000001

void presentacion()
{
   printf("buenos dias\n");	
}

bool comprobadorError(int cant_errores, float errores[]) {
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
    int cant_var;
    printf("\nCantidad de variables del sistema: ");
    scanf("%d", &cant_var);
    int tam_matriz = cant_var; // No necesitas aumentar en 1 aquí

    float matriz[tam_matriz][tam_matriz + 1]; // Columna adicional para el término independiente

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

    // Método de Jacobi
    float errores[cant_var];
    float x[cant_var];
    float prev_x[cant_var];
    for (i = 0; i < cant_var; i++) {
        errores[i] = 1; // Inicializa errores grandes
        x[i] = 0; // Valores iniciales
    }

    int contador = 0;
    printf("Iteración | ");
    for (i = 0; i < cant_var; i++) {
        printf("x[%d]       | ", i + 1);
    }
    printf("Errores\n");
    printf("----------------------------------------------------------\n");

    while (!comprobadorError(cant_var, errores)) {
        // Impresión como en Excel
        printf("%-10d | ", contador);
        for (i = 0; i < cant_var; i++) {
            printf("%-10f | ", x[i]);
        }

        // Actualización de los valores
        for (i = 0; i < cant_var; i++) {
            prev_x[i] = x[i]; // Guardar valor anterior
            float suma = 0;

            // Sumar los términos de la ecuación
            for (j = 0; j < cant_var; j++) {
                if (j != i) {
                    suma += matriz[i][j] * prev_x[j]; // Usa el valor anterior
                }
            }

            // Calcular el nuevo valor
            x[i] = (matriz[i][tam_matriz] - suma) / matriz[i][i]; // Término independiente

            // Calcular el error
            errores[i] = fabs(x[i] - prev_x[i]);
            prev_x[i] = x[i]; //profe si está leyendo hasta acá está linea no salía en google y casi
                           //pierdo mi estabilidad mental :'(
        }

        // Mostrar errores
        for (i = 0; i < cant_var; i++) {
            printf("%-10f | ", errores[i]);
        }
        printf("\n");

        contador++; // Incrementar el contador de iteraciones
    }

    // Mostrar resultados finales
    printf("\nResultados finales:\n");
    for (i = 0; i < cant_var; i++) {
        printf("x[%d] = %f\n", i + 1, x[i]);
    }
}

void doolittle() {
    int i, j, k;
    int cant_var;
    printf("\nCantidad de variables del sistema: ");
    scanf("%d", &cant_var);
    int tam_matriz = cant_var;

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
void newtonRaphson()
{
	printf("Trabajando en el metodo de newton Raphson\n");
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