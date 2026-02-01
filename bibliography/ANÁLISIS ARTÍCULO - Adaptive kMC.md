# ANÁLISIS ARTÍCULO - Adaptive kMC model



## 1. Introducción



* Se desarrolla un modelo microscópico que utiliza una tasa de adsorción especializada que se adapta a diferentes regímenes de crecimiento considerando los efectos termodinámicos de la superficie para diferentes mecanismos de crecimiento. Específicamente, la tasa de adsorción adaptativa redefine la fuerza impulsora de la cristalización al considerar las energías de adhesión de diferentes sitios de adsorción y la sobresaturación.

* En el caso de cristales creciendo en solución, el grado de sobresaturación emerge como la variable principal que controla el mecanismo de crecimiento. En general, existen tres regímenes importantes: espiral, escalonado y rugoso.

* Primero, en el caso de un régimen de crecimiento en espiral a baja sobresaturación, prevalece el crecimiento de defectos, impulsado por la mayor probabilidad de adherencia debido a la presencia de defectos superficiales. Estos defectos, reducen la barrera de energía de adhesión, acelerando así el crecimiento del cristal a medida que avanza el escalón. Con una deposición adicional, el crecimiento del cristal continúa y emerge un nuevo escalón, que crece a lo largo de una nueva trayectoria, pero aún paralelo a la superficie. A medida que este segundo paso se propaga, deja el defecto a su paso, lo que da lugar a otro paso. A medida que esto continúa, el efecto neto es el crecimiento de cristales que giran en espiral alrededor del defecto del tornillo.

* Segundo, en el régimen de crecimiento por pasos (o escalonado) a mayores niveles de sobresaturación, la deposición ocasional sobre superficies cristalinas planas se hace posible a medida que el aumento de la sobresaturación disminuye la energía libre hasta un punto en el que la deposición puede ocurrir incluso en una superficie limpia y plana. Esto conduce a un crecimiento bidimensional, también conocido como crecimiento de "nacimiento y propagación", ya que la deposición posterior tiende a ocurrir en los pasos recién formados.

* En tercer lugar, a niveles de sobresaturación relativamente altos, el crecimiento cristalino se produce en un régimen de crecimiento irregular (rugoso). A niveles de sobresaturación más altos durante el crecimiento cristalino, las unidades de crecimiento tienden a adsorberse en adátomos debido a la mayor disponibilidad de sitios de adsorción y a la rapidez del proceso de crecimiento. Estos adátomos no contribuyen al crecimiento de superficies cristalinas bien estructuradas, sino que conducen a la formación de islas aisladas. Con el tiempo, estas islas pueden fusionarse y coalescer a medida que evolucionan las condiciones de crecimiento, desarrollando finalmente una capa cristalina continua.

* Además de los mecanismos de crecimiento descritos dentro de los cristalizadores, es esencial tener en cuenta una característica específica conocida como el hang en la evolución de la morfología superficial de los cristales. El hang es un sitio de adsorción distintivo que se caracteriza por la presencia de tres vecinos más cercanos existentes. La formación del hang típicamente implica la alineación de un pliegue positivo o negativo y un escalón en la superficie del cristal. Esta configuración crea un sitio con mayor potencial de adsorción debido a la disponibilidad de estructuras vecinas, lo que promueve la adherencia de las unidades de crecimiento. El hang desempeña un papel crucial en la dinámica del crecimiento cristalino, proporcionando una faceta adicional a la compleja interacción de las estructuras superficiales y los niveles de sobresaturación. Representa un escenario de adsorción único que contribuye a las complejidades generales de los mecanismos de crecimiento cristalino, influyendo aún más en la evolución de la morfología a medida que varían las condiciones en el cristalizador.

* Las simulaciones kMC presentan una vía alternativa debido a su capacidad para abarcar dimensiones de red más grandes (por ejemplo, y escalas de tiempo extendidas, por ejemplo, que abarcan desde 1 µs hasta 10 000 s).



---



## 2. Descripción microscópica del crecimiento cristalino



* Los regímenes de crecimiento cristalino se refieren a diferentes mecanismos y patrones de crecimiento observados durante la formación de cristales. Estos regímenes están influenciados por la interacción entre las unidades de crecimiento del cristal y la superficie sobre la que crecen. De ahí los tres principales regímenes de crecimiento.

* La cristalización en solución es el proceso mediante el cual las unidades de crecimiento parten de un estado líquido desordenado de mayor energía y transicionan a un estado sólido menos energético y ordenado. Esta transición de fase es dirigida por un decrecimiento del potencial químico. La tasa a la que esta transición ocurre se conoce como tasa de crecimiento del cristal, y depende de la diferencia del potencial químico entre los estados desordenados y ordenados.

* La energía requerida para que las unidades de crecimiento se unan a la superficie del cristal, conocida como barrera de energía libre, también contribuye a esta diferencia de potencial químico. Esta barrera depende de la densidad de los sitios de adsorción, a lo largo de las caras del cristal donde pueden unirse nuevas unidades de crecimiento.

* A mayor sobresaturación, la diferencia de potencial químico es mayor, lo que proporciona una mayor fuerza impulsora para la cristalización. Esto conduce a mayores tasas de crecimiento de los cristales. Cuando la tasa de crecimiento es muy alta en condiciones de alta sobresaturación, los cristales tienden a crecer de forma irregular e incontrolada en lugar de formar superficies lisas y estructuradas. Este crecimiento irregular (rugoso) se produce porque los UC se unen aleatoriamente a los adátomos de la superficie en lugar de hacerlo preferentemente a escalones o pliegues.

* Los adátomos actúan como sitios de nucleación para la formación de islas cristalinas aisladas. Estas islas pueden unirse y fusionarse con el tiempo, formando finalmente una capa cristalina continua si las condiciones de crecimiento cambian o se estabilizan. Generalmente se observa que cuando la tasa de crecimiento de los cristales es alta, esto conduce a la nucleación y al crecimiento de estas islas de cristales individuales que no forman una capa de cristal continua.

* La tasa de crecimiento de los cristales depende en gran medida de la sobresaturación, ya que la principal fuerza impulsora de la cristalización es el potencial químico.

* De esta forma, se observa un crecimiento rugoso para valores altos de sobresaturación durante la cristalización en solución.

* En niveles de sobresaturación más bajos, el crecimiento de nucleación bidimensional se convierte en el mecanismo dominante para el crecimiento de los cristales.

* La nucleación bidimensional requiere una barrera energética menor que el crecimiento rugoso, lo que le permite predominar en sobresaturaciones más bajas. Las nuevas UC se unen a los bordes de las capas cristalinas existentes, formando terrazas o mesetas planas en la superficie del cristal. Procede mediante la formación de pequeños núcleos en la superficie, que luego se expanden lateralmente. Esto resulta en un crecimiento capa por capa y una morfología plana. La superficie cristalina se caracteriza por escalones claramente definidos que representan los límites entre las diferentes capas cristalinas. Este tipo de crecimiento se relaciona con sitios de adsorción bordes (edges) con un solo vecino. Cuando las UC se adsorben en estos sitios de borde, se incorporan a la red cristalina y provocan el crecimiento de nuevas capas paralelas a la superficie cristalina existente. El crecimiento se produce mediante la adición de capas, una sobre otra, lo que resulta en la formación de terrazas de espesor uniforme.

* Con una sobresaturación muy baja, el crecimiento de los cristales se vuelve insignificante y solo puede continuar cuando la energía de adhesión del sitio es alta. Los sitios kink son muy importantes porque las UC que se unen allí forman más enlaces con las UC vecinos que las que se unen a las terrazas o a los bordes planos de los escalones. Estos kinks (puntos de quiebre) son relativamente estables y actúan como puntos de nucleación para un mayor crecimiento cristalino. Los kink son lugares donde termina un escalón cristalino, con dos UC vecinas.

* A medida que las UC se adsorben en los kink, contribuyen al crecimiento del cristal en espiral. El crecimiento en espiral se caracteriza por la adición continua de nuevas UC al cristal en una disposición helicoidal o espiral. Las UC recién añadidas se unen a los kink y se propagan a lo largo de la superficie del cristal, formando una espiral continua. El crecimiento en espiral es especialmente frecuente cuando el cristal está expuesto a bajos niveles de sobresaturación o en presencia de defectos superficiales específicos, como la dislocación helicoidal.

* La interacción entre las UC y los sitios microscópicos influye en el comportamiento de crecimiento, la morfología del cristal y las propiedades del cristal resultante.

* La energía de unión se refiere a la energía de enlace que se libera cuando una UC se une a la superficie de una cara cristalina.

* Hay indicios de una relación entre la energía de unión y la velocidad a la que una UC se une a una cara cristalina. Se sugiere que la energía de enlace es inversamente proporcional al tiempo necesario para formar un enlace. A medida que el tiempo necesario para la formación del enlace disminuye con el aumento de la energía de unión, la velocidad de desplazamiento de una cara cristalina también aumenta. Esto implica que la cara cristalina con mayor energía de unión crece más rápido, lo que conduce al equilibrio de la forma cristalina durante el crecimiento cristalino.

* Por lo tanto, la tasa de crecimiento del cristal depende de la energía de adhesión de los GU del cristal. La energía de unión se convierte en un parámetro representativo de importancia morfológica, ya que las caras cristalinas con mayor energía de unión dominan el régimen de crecimiento y juegan un papel crucial en la conformación de la morfología superficial final del cristal.

* El modelo kMC desarrollado se puede adaptar a diversas estructuras cristalinas y materiales en función de sus energías de unión (adhesión) para su respectiva UC.



---



## 3. Modelo kMC desarrollado



* Basado en la teoría de transiciones de estado (TST) para simular dinámicas discretas de largo plazo de forma eficiente. Cada estado corresponde a un mínimo energético local. A cada posible transición se le asigna una probabilidad (constante de velocidad o de tasa) basada en la TST. La TST proporciona una forma aproximada de calcular la constante de velocidad para cada posible transición entre estados y permite la utilización de factores externos e internos que afectan la cristalización. La TST supone que existe una superficie divisoria que separa los estados y calcula la velocidad como el flujo de equilibrio a través de esta superficie divisoria. Esto evita tener que simular explícitamente las transiciones.

* Este modelo ignora las diferencias en la frecuencia de vibración y la energía libre rotacional de las moléculas en diferentes posiciones sobre la superficie del cristal.

* En el estudio del artículo, la red cristalina se inicializa aleatoriamente para tener una mezcla de diferentes configuraciones de GU para eliminar cualquier sesgo hacia un solo tipo de configuración. Además, se ignoran la convección, el calor latente de cristalización y las interacciones laterales.

* En la cinética superficial, la fuerza impulsora detrás de este proceso de cristalización a menudo se ha atribuido al cambio general en la energía libre superficial durante la cristalización, junto con el cambio en el potencial químico de la solución que experimenta la cristalización. Cuanto mayor sea la magnitud del diferencial del potencial químico, más pronunciado será el crecimiento del proceso de cristalización. (Eq 1.)

* Con la adición de unidades de crecimiento a una superficie de cristal plana macroscópica, hay un cambio en la energía libre de Gibbs por mol (a T y P fijas) de todo el sistema correspondiente a (Eq 2.)

* El tamaño crítico correlaciona la morfología de la superficie con la supersaturación y predice los diferentes regímenes de crecimiento observados. Cuando la sobresaturación está en un régimen de crecimiento en espiral, el tamaño crítico se hace grande y observamos una superficie más lisa a menos que haya un defecto superficial en la red. A medida que la sobresaturación aumenta al régimen de crecimiento escalonado, se crea un grupo de islas con un tamaño aproximadamente igual a su tamaño crítico. A medida que la sobresaturación aumenta al régimen de crecimiento rugoso, el tamaño de la isla disminuye hasta que es más pequeño que el tamaño de la unidad GU, y comienza la rugosidad de la morfología de la superficie. Si un grupo es más grande que este tamaño crítico, tiene una mayor probabilidad de crecer que de decaer, mientras que si el grupo es más pequeño que este tamaño crítico, tiene una mayor probabilidad de decaer. (Eqs 3 y 4.)

* Para integrar este mecanismo en el modelo kMC, primero necesitamos calcular el cambio de energía libre para cada UC transferido durante el evento de adsorción, dado por (Eq 5.)

* Consideraciones de la (Eq 5.) hacen postular que el tiempo requerido para la formación de un enlace en un sitio de adsorción es inversamente proporcional a la energía del enlace. Como resultado, la velocidad de desplazamiento de una cara de cristal aumenta con una mayor energía de fijación (unión, adhesión), lo que hace que la energía de unión sea un indicador significativo de su importancia morfológica. En otras palabras, una cara de cristal desaparecerá antes durante el crecimiento del cristal si se aleja del centro de crecimiento a un ritmo más rápido.

* Lo anterior indica que la tasa de adsorción de cada UC depende de la barrera de energía libre de la siguiente manera: (Eqs 6 y 7.)

* La superficie del cristal se encuentra en equilibrio dinámico con las UC en la solución circundante. Las UC se adhieren constantemente a la superficie del cristal (adsorción) y se desprenden de la superficie de nuevo a la solución (desorción). El crecimiento a partir de una solución sobresaturada se produce porque el flujo de UC adheridas a la superficie del cristal supera el flujo de UC desprendidas de la superficie.

  

  ---

  

  ## 4. Resultados y discusión

  

* 
