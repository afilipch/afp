{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Microarray analyses for C.glutamicum"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "%matplotlib inline\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "We have already processed the initial table in a way that: \n",
    "* it is filled by numeric values\n",
    "* formated to tab-separated version (.tsv)\n",
    "\n",
    "Now we download the table as DataFrame:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "PATH_TABLE = '/hdd/projects/angela_metaanalyses/data/preprocessed.tsv'\n",
    "df = pd.read_csv(PATH_TABLE, index_col = 0, sep='\\t')\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    " ## Basic statistics :"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Number of experiments:\t403\n",
      "Number of genes:\t3047\n"
     ]
    }
   ],
   "source": [
    "print(\"Number of experiments:\\t%d\\nNumber of genes:\\t%d\" % df.shape)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Mean microarray intensity distribution along the Experiments:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 864x288 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "def draw_mean_intensity_histogram(df, axis_num, percentile, label, bins):\n",
    "    exp_means = df.mean(axis=axis_num)\n",
    "    resolution_limit = np.percentile(exp_means, percentile)\n",
    "    exp_means_resolution = [x for x in exp_means if x <= resolution_limit]\n",
    "\n",
    "    fig, axes = plt.subplots(1,2, figsize=(12,4))\n",
    "    fig.suptitle(label)\n",
    "    for ax, values, title in zip(axes, [exp_means, exp_means_resolution], ['all', 'top %d' % percentile]):\n",
    "        ax.set_ylabel(\"Number of %s\" % label)\n",
    "        ax.set_xlabel(\"Mean Expression\")\n",
    "        ax.hist(values, bins = bins);\n",
    "\n",
    "draw_mean_intensity_histogram(df, 1, 90, \"Experiments\", 25)\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Mean microarray intensity distribution along the Genes:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 864x288 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "draw_mean_intensity_histogram(df, 0, 90, \"Genes\", 100)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "We split the DataFrame into:\n",
    "* microarray intensities (DataFrame values) as numpy array\n",
    "* experiment names (DataFrame rows) as list\n",
    "* gene names (DataFrame columns) as list"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
