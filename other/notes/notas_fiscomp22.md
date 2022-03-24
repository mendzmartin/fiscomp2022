# Física Computacional
## 15/03/2022

Aula 32 -> Práctico
Aula 13 -> Teórico

### Contenido que vamos a ver
+ Linux
+ Shell, es decir, un software que controla la terminal. En particular veremos `bash-shell`.
+ Graficador. En particular, veremos `gnuplot` ya que es un graficador liviano y útil debuggear código.
+ Lenguaje de programación. En particular, usaremos `FORTRAN`, y aprenderemos el lenguaje indirectamente, resolviendo problemas físicos concretos. Sin embargo, se requiere cierto conocimiento previo del lenguaje.

### Linux

El comando `ls -a` muestra en formato tipo lista los archivos dentro de un directorio.
El comando `ls -al` muestra en formato tipo lista los archivos dentro de un directorio y detalles.
El comando `ls -alh` muestra en formato tipo lista los archivos dentro de un directorio. Aquí `h` significa human, cambia las unidades para que sean comodas)

Con el comando `ls -alh` desde la terminal tendremos una salida de la forma,
```
	drwxrwxr-x 2 mendez mendez 4,0K mar 15 09:42 .
```
donde:
+ 1er columna `drwxrwxr-x` nos muestra los permisos que tienen los usuarios (dueño, grupo o invitados).
+ 2da y 3er columna nos muestra el dueño y el grupo al que pertenece el archivo.
+ 4ta columna nos muestra el tamaño del archivo.
+ 5ta, 6ta y 7ma columna nos muestra el dia, la fecha y la hora de la última modificación del archivo.

Respecto a los permisos notemos que:
+ letra `d` -> permiso para crear
+ letra `r` -> permiso para leer
+ letra `w` -> permiso para escribir
+ letra `x` -> permiso para ejecutar
+ simbolo `-` -> no se tiene permiso

Además:
+ Los primeros tres caracteres son los permisos que tiene el dueño
+ Los segundos tres caracteres son los permisos que tiene el grupo
+ Los terceros tres caracteres son los permisos que tiene el visitante

Los archivos que empiezan con un punto `.` son archivos ocultos de linux.

#### Comandos varios

+ `alias` nos muestra los alias definidos en el bash shell.
+ `ls P*` nos muestra una lista de todos los directorios que empiezan con el caracter `P`. Hay que tener cuidado porque el símbolo `*` es un *comodín* y no usar nunca el comando `rm *` para remover algun archivo porque borra todo y sin posibilidad de recuperación (o al menos muy dificil de recuperar).
+ `rm -i archivo` lo podemos usar para borrar archivos y nos va a preguntar si estamos seguros de querer borrarlo. Es más, podemos ejecutar el comando `alias rm = "rm -i"` para que cada vez que corramos el comando `rm` nos ejecute el comando alias `rm -i`. Sin embargo, esto se debe hacer cada vez que abrimos una terminal lo cual no es muy práctico, para que sea automático debemos modificar el archivo `.bashrc`.
+ `chmod {u,g,o}{+-}{dwrx} archivo.extension` para cambiar los permisos de un archivo donde:
  +  `+` significa agregar permiso y `-` significa quitar permiso.
  + `u` (user) se refiere al usuario o dueño, `g` se refiere al grupo (group) y `o` se refiere a otros (other)
+ `jobs` para ver los procesos que se estan ejecutando en el background.
+ `fg` para remover un proceso del background y traerlo a la terminal.

### visor de texto
#### more
+ Para ejecutar corremos el comando `more mynote.txt` y dentro de la visualización podemos escribir `/palabra` para que nos muestre la palabra buscada en el documento y con la tecla `n` podemos ir mostrando todos los resultados encontrados.

### Editores de Texto
#### Kate
+ Para poder abrir Kate desde la terminal y permitir seguir trabajando en la terminal, es decir, pedir que Kate se ejecute en el background tendremos que hacer
```
	$ kate 2> /dev/null 
	$ [Ctrl+z]
	[1]+  Detenido                kate 2> /dev/null
	$ bg
	[1]+ kate 2> /dev/null &
```
#### nano

#### evince
+ Para ejecutar hacemos
```
$ evince archivo.pdf &
```
