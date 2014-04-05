# $1 := type (for now, either state or control)
# $2 := exec file
# $3 := timestep
# $4 := symeval c dir
# $5 := masks dir
# $6 := sym func params

TYPE=$1
EXEC=$2
T=$3
SYM_DIR=$4
MASK_DIR=$5
SYM_FUNCS=$6

DES_C="${SYM_DIR}/${TYPE}-symeval-${T}.c"
DES_TXT="${MASK_DIR}/${TYPE}-masks-${T}.txt"

TMP=tmp
TMP_C_BODY=tmp_c_body
echo " " > $TMP_C_BODY # to initialize

# add new file top
SRC_C_TOP=$'#include "'$TYPE'-symeval.h"
#include <math.h>
double Square(double x) {
    return x*x;
}'

# move to correct folder
echo "$SRC_C_TOP" > $TMP

for SYM_FUNC in $SYM_FUNCS
do
    echo ${SYM_FUNC}

    rm -f *.c, *.txt

    mono $EXEC $T $SYM_FUNC

    SRC_C=$(ls *.c | head -1)
    SRC_TXT=$(ls *.txt | head -1)

    # hard code delete top
    sed -i '1,71d' $SRC_C
    # delete from first }; to the end
    START_TAIL_CUT=$(grep -n "};" $SRC_C | cut -f1 -d: | head -1)
    sed -i ${START_TAIL_CUT}',$d' $SRC_C
    # replace function definition name
    sed -i "s/eval/eval${SYM_FUNC}/" $SRC_C

    cat $SRC_C >> $TMP_C_BODY

    # copy to correct folder
    cp $SRC_TXT $DES_TXT

done

# prepare masks for insertion into top of c file
NUM_VARS=$(wc -w $SRC_TXT | cut -f 1 -d " ")
#sed -i "s/ /, /g" $SRC_TXT
MASK=$(sed -n '1p' $SRC_TXT)

SRC_C_MASK=$'
char* getMask() {
    return "'$MASK'";
}'

#    int mask['$NUM_VARS'] = {'$MASK'};

echo "$SRC_C_MASK" >> $TMP

# add c src to rest of tmp
cat $TMP_C_BODY >> $TMP

cp $TMP $DES_C
rm -f $TMP, $TMP_C_BODY, *.c, $EXEC, *.txt
