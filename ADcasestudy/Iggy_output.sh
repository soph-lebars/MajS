
# Will give the output file of iggy (predicted node) for the observed file passed in argument 

var1=$1

~/Documents/These/iggy/target/release/iggy -n regfromdag.cif -o $var1 -a -p |   sed -n '/Predictions/, /Prediction/{ /Predictions/! { /Prediction/! p } }'  | sed -r '/^\s*$/d' > "Iggy_${var1/obs/txt}"

