# Análisis del Proyecto

## Objetivo Principal

El objetivo de este proyecto es implementar un sistema de **renderizado no fotorrealista (NPR)**. Concretamente, se busca transformar una imagen de entrada (por ejemplo, una pintura famosa) en una versión estilizada compuesta por una serie de pinceladas (`Strokes`).

El código base proporcionado actúa como un motor de renderizado, manejando las operaciones de bajo nivel como la carga de imágenes, la creación de un lienzo (`Canvas`) y el dibujo de pinceladas individuales. La tarea principal del desarrollador es **implementar la lógica algorítmica** para determinar el conjunto óptimo de pinceladas que mejor se aproxime a la imagen objetivo. Esto convierte el problema en un desafío de **optimización**.

## Componentes Clave

- **Lienzo (Canvas):** Representa la imagen RGB que se está pintando. Se inicializa en blanco y se va modificando con cada pincelada.

- **Pinceles (Brushes):** Son imágenes en escala de grises que funcionan como máscaras para dar textura y forma a las pinceladas.

- **Pincelada (Stroke):** Es la unidad fundamental de dibujo. Cada `Stroke` es una estructura de datos con atributos que definen su apariencia y posición:
    - Posición (`x_rel`, `y_rel`)
    - Tamaño (`size_rel`)
    - Rotación (`rotation_deg`)
    - Color (`r`, `g`, `b`)
    - Tipo de pincel (`type`)

- **Motor de Renderizado:** Las funciones en `stroke.cpp` y `stroke.h` permiten dibujar un `Stroke` en el `Canvas`, aplicando transformaciones (escala, posición, rotación) y mezclando el color según la máscara del pincel.

## Flujo de Ejecución Típico

1.  **Carga:** Se carga una imagen objetivo (ej. `instancias/mona.png`) en un `Canvas` para usarla como referencia.
2.  **Generación de Solución:** Se crea un conjunto (un `std::vector`) de `Stroke`. Aquí es donde reside el núcleo del algoritmo a implementar. El objetivo es generar los `Strokes` que recreen la imagen objetivo.
3.  **Renderizado:** Se renderiza el vector de `Strokes` sobre un `Canvas` en blanco para crear la imagen final.
4.  **Evaluación:** Se compara el `Canvas` generado con el `Canvas` de la imagen objetivo. El `README` sugiere usar el **Error Cuadrático Medio (MSE)** como métrica para cuantificar la similitud.
5.  **Optimización:** Se modifican los atributos de los `Strokes` (su posición, color, tamaño, etc.) de forma iterativa para minimizar el error (MSE) y así mejorar la calidad de la pintura generada.

En resumen, el proyecto es un esqueleto para un programa de pintura automática donde el desafío principal es diseñar e implementar un algoritmo de optimización para encontrar la mejor representación de una imagen usando un número limitado de pinceladas.
