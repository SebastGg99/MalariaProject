# notebooks - Guia de experimentos reproducibles

Este directorio contiene notebooks de exploracion, calibracion y reproduccion de resultados del modelo kMC.

## Inventario de notebooks

- `Malaria002D_v1.ipynb`: experimentos del modelo base 2D.
- `paper_kMC.ipynb`: reproduccion y analisis de la variante del paper.
- `experiments_1.ipynb` a `experiments_4.ipynb`: baterias de pruebas y sensibilidad de parametros.

## Archivos auxiliares

- `beta-hematina.txt`: datos de referencia usados en algunas comparaciones.

## Orden sugerido de ejecucion

1. Cargar entorno y dependencias del proyecto.
2. Ejecutar celdas de import y configuracion de ruta.
3. Ejecutar bloques de inicializacion de parametros y red.
4. Correr simulaciones por bloques de tiempo.
5. Generar graficas de crecimiento y probabilidad de sitio.

## Buenas practicas

- Fijar semillas cuando se requiera reproducibilidad.
- Evitar mezclar resultados de experimentos con parametros distintos en una misma sesion.
- Exportar figuras finales a `results/` con nombres descriptivos.

## Notas para Sphinx

Este README puede incluirse en la documentacion de uso para explicar el flujo experimental sin duplicar contenido.
