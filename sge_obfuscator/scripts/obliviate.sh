#!/usr/bin/env bash

#
# Change these as needed
#
SAGE='sage'
PYTHON='python2'
CODE_DIR='../code'
CIRCUIT_DIR='../code/circuits'
LOG_DIR='runs'
SCHEME=''

mkdir -p $LOG_DIR

secparam=24

circuit="fourandsnot.circ"
eval='00000'
dir="$LOG_DIR/obliviate.$circuit.$secparam"
mkdir -p $dir


echo "* Running $circuit with security parameter $secparam (oblivious)"
$SAGE $CODE_DIR/obfuscator obf \
    --load-circuit $CIRCUIT_DIR/$circuit \
    --secparam $secparam \
    --obliviate \
    $SCHEME \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obv-obf-time.log
obf=$circuit.obf.$secparam
du --bytes $CIRCUIT_DIR/$obf/* \
    | tee $dir/$circuit-$secparam-obv-obf-size.log
$SAGE $CODE_DIR/obfuscator obf \
    --load-obf $CIRCUIT_DIR/$obf \
    --eval $eval \
    $SCHEME \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obv-eval-time.log
rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam

echo "* Running $circuit with security parameter $secparam (std)"
$SAGE $CODE_DIR/obfuscator obf \
    --load-circuit $CIRCUIT_DIR/$circuit \
    --secparam $secparam \
    $SCHEME \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obf-time.log
obf=$circuit.obf.$secparam
du --bytes $CIRCUIT_DIR/$obf/* \
    | tee $dir/$circuit-$secparam-obf-size.log
$SAGE $CODE_DIR/obfuscator obf \
    --load-obf $CIRCUIT_DIR/$obf \
    --eval $eval \
    $SCHEME \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-eval-time.log
rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam

circuit='fourxors.circ'
echo "* Running $circuit with security parameter $secparam (oblivious)"
$SAGE $CODE_DIR/obfuscator obf \
    --load-circuit $CIRCUIT_DIR/$circuit \
    --secparam $secparam \
    --obliviate \
    $SCHEME \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obv-obf-time.log
obf=$circuit.obf.$secparam
du --bytes $CIRCUIT_DIR/$obf/* \
    | tee $dir/$circuit-$secparam-obv-obf-size.log
$SAGE $CODE_DIR/obfuscator obf \
    --load-obf $CIRCUIT_DIR/$obf \
    --eval $eval \
    $SCHEME \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obv-eval-time.log
rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
