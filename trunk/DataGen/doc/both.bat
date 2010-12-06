datagen geodb 400 1 5 | rotateq .998 0 -.063 0   | rotateq .998 0 -.063 0   | rotateq .998 -.063 0 0   | transl .8 -.5 -1 > geodb.csv
datagen ranusp 400  | scale .5  > sphere.csv
copy/b geodb.csv+sphere.csv both.csv
type both.csv | rotateq .998 0 -.063 0   | rotateq .998 0 -.063 0  |svg 100 > bothROT.svg
type both.csv | rotateq .998 0 -.063 0   | svg 100 > both.svg
