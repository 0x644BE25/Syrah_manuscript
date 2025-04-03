# Syrah manuscript
Post-Syrah data and analysis/plot generation code for the manuscript [*Syrah: an analysis pipeline to maximize spatial transcriptomic output*](URLRULRUL)
The [main Syrah pipeline GitHub repostiory is here](https://github.com/0x644be25/Syrah) and contains a [Tutorial](https://github.com/0x644BE25/Syrah/blob/main/Syrah_tutorial.md) to walk you through installing and running Syrah.

# Generating plots
The R script [`generate_plots_and_xlsx.R`]() uses the [`data`](https://github.com/0x644BE25/Syrah_manuscript/tree/main/data) folder and contained files as input and generates all plots and underlying data from the maunuscript. If you have run Syrah on the raw data yourself (see below) you must ensure that the filenames match those in the [`data` folder](https://github.com/0x644BE25/Syrah_manuscript/tree/main/data) but you **do not** need to filter based on # UMIs, as the plotting script will do that for you.

# Running Syrah
If you would like to run Syrah on the raw data yourself, these instructions for [running Syrah](https://github.com/0x644BE25/Syrah_manuscript/blob/main/running_Syrah.md) will guide you through installation and processing all three datasets from the manuscript.

# Troubleshooting
The main Syrah repo has more information, but if you're still having problems don't hesitate to contact me through [email](cbrewster@stowers.org) or on GitHub. The whole point of making Syrah is for *YOU* to be able to use it!
