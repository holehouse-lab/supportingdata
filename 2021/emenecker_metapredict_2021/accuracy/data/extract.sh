#!/bin/bash

echo "Extracting data..."

tar -xvf data1.tgz
mv data1/* .

tar -xvf data2.tgz
mv data2/* .

rmdir data1
rmdir data2

echo "Data extracted!"

