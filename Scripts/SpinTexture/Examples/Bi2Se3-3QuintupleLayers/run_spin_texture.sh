# Run all calculations
for x in 'bands' 'circle' 'grid' 'constant_energy'; do
    echo "Calculation with mode ${x}"
    python ../../SpinTexture.py "input.${x}.json"
done

# Plot all results
for x in 'bands' 'spinbands' 'circle' 'grid' 'constant_energy'; do
    echo "Plotting: ${x}"
    python ../../PlotSpinTexture.py "plot.${x}.json"
done
