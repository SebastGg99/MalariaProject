# data - Convenciones de datos

Este directorio esta reservado para datos de entrada del proyecto.

## Estado actual

Actualmente no hay datasets versionados en esta carpeta.

## Que tipos de datos van aqui

- Tablas experimentales de referencia.
- Series temporales para validacion.
- Insumos para calibracion de parametros.

## Reglas de formato recomendadas

- Preferir CSV o TXT tabular con encabezados claros.
- Incluir unidades en nombres de columnas cuando aplique.
- Evitar archivos binarios sin documentacion asociada.

## Trazabilidad

Para cada dataset nuevo, documentar:

- Fuente.
- Fecha.
- Preprocesamiento aplicado.
- Script o notebook que lo consume.

## Notas para Sphinx

Este README debe enlazarse desde la seccion de uso o metodologia para mantener trazabilidad de insumos.
