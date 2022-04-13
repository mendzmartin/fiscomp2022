#!/bin/bash

#################################################################
### Purpose: script para correr en cluster bandurria
#################################################################
#################################################################
### COMANDOS ESPECÍFICOS PARA ADMINISTRADOR DE COLAS SGE
#################################################################

# Current Working Directory. Directorio actual es el raiz.
# Se trabajará por defecto en el directorio indicado,
#   con los ficheros incluidos en él.
#$ -cwd

# Manda un email si pasa algo con el proceso
##$ -m eas -M martinmendez@mi.unc.edu.ar
#$ -M martinmendez@mi.unc.edu.ar 
#$ -m bea

# Nombre del proceso.
# Por defecto toma el nombre del script que se ejecuta.
#$ -N heylchaos

# Especifica si la salida de error de la tarea a ejecutar va a ser
#   la salida de error de SGE, es decir, si stdout y stderr apuntan
#   al mismo archivo de salida.
#$ -j y

# Usar shell específico (en este caso bash) para el proceso a correr
#$ -S /bin/bash

# Pido memoría RAM para el proceso
#  si quiero 4GB reservo un poco más
#  (por ej. si necesito 2GB recervo 2.2GB)
#$ -l mem_free=4.2G

# Obligatorio en procesos paralelos (-pe [entorno_paralelo] [cpus])
#  donde [entorno_paralelo] = [smp/make/openmpi/lam/mpich]
#  donde [cpus] = número de cpus
#$ -pe smp 1

# NSLOTS -> Variable de entorno que especifica el número de slots
#           de la cola asignados al trabajo paralelo.
#$ -v OMP_NUM_THREADS=${NSLOTS}

# Reservar slots a medida que otros procesos los liberan
#$ -R y

# Tiempo de CPU (wall clock) que se solicita para el proceso
#  h_rt=[h:min:sec]
#$ -l h_rt=4:00:00

# Especificar que el proceso es capaz de hacer checkpoints
#$ -ckpt dmtcp

# Eligir el archivo output.log para que allí se escriba
#  el standard output (salida)
#$ -o output_run.log

# Excluye este/estos nodos
#

#################################################################
### EXPORTACIÓN DE VARIABLES DE ENTORNO
#################################################################

./script_run.sh

#################################################################
### NOTAS Y/O COMENTARIOS
#################################################################

# " #$ " simple numeral con pesos es un comando
# " #  " numeral solo es un comentario
# dentro de bandurria usar el comando [/home/user]
# 	qsub sge_script.sh 	-> para encolar el trabajo
# 	qstat 			-> para ver los trabajos encolados
# 	qdel jobID		-> para remover un trabajo de la cola
#	qacct -j jobID 		-> report and account for Sun Grid Engine usage
#	qdel -u username
