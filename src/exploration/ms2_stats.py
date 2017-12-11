"""
- how many MS2's do we have per molecule?
    - how many peaks per MS2 (per molecule)?
    - basically: how much of our data is "bad"?
    - is this biased across classes?
-how many molecules of each class are in each type of ionization mode?
-Binning

Initialize ionization type dictionary {class: [positive_count, negative_count, N/A count, None count],...}
Initialize MS2/molecule count list {MS2/molecule#, ...}
Initialize peaks/ms2 count list {peaks/MS2#,...}
Initialize [(# of MS2's/taxonomy), peaks/MS2] dictionary: {superclass: MS2/molecule#, kingdom: MS2/molecule#, ...}

For each metabolite:
  update ionization dictionary/taxonomy level: ['Positive', 'Negative', 'N/A', None]
  append to MS2/molecule count list
  append to peaks/MS2 count list
  [# of MS2's/taxonomy level, running avg peaks/MS2, count] dictionary at each taxonomy level


Generate histogram of MS2/molecule
Generate histogram of peaks/MS2
Generate 4 bar graphs of average peaks/MS2 for each [class, subclass, kingdom, superclass]
Generate 4 bar graphs of average MS2's/molecule for each [class, subclass, kingdom, superclass]



Find min/max m/z
Generate bins

"""

import os
import sys
import csv
sys.path.append(os.getcwd() + '/../util/')
import util
import pdb
import matplotlib.pyplot as plt
import sys
import numpy as np

sys.path.insert(0, '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
import pandas as pd
import seaborn as sns


class StatsMaker:
  def __init__(self):
    csv_path = os.getcwd() + '/../../data/clean/metabolites_and_spectra.csv'
    self.metabolite_dict = util.unpackCSV(csv_path)
    self.kingdom_ionization_dict = dict()
    self.sub_class_ionization_dict = dict()
    self.super_class_ionization_dict = dict()
    self.class_ionization_dict = dict()
    self.ms2_per_molecule = list()
    self.peaks_per_ms2 = list()
    self.ionization_modes = ('Positive', 'Negative', 'N/A', None)
    self.kingdom_dict = dict()
    self.sub_class_dict = dict()
    self.super_class_dict = dict()
    self.class_dict = dict()


  def updateIonization(self, taxonomy_dict, ms2_list):
    if not taxonomy_dict:
      return
    for ms2 in ms2_list:
      ionization_mode = ms2.ionization_mode
      index = self.ionization_modes.index(ionization_mode)

      for (dictionary, tax) in [(self.kingdom_ionization_dict, 'kingdom'), (self.super_class_ionization_dict, 'super_class'), (self.sub_class_ionization_dict, 'sub_class'), (self.class_ionization_dict, 'class')]:
        tax_instance = taxonomy_dict.get(tax)
        if tax_instance:
          if not dictionary.get(tax_instance):
            dictionary[tax_instance] = [0, 0, 0, 0]
          dictionary.get(tax_instance)[index] += 1


  #dictionaries at each taxonomy level have format: {kingdom_name: [total # of MS2's, running avg, count]}
  def updateCountingDicts(self, taxonomy_dict, ms2_list):
    if not taxonomy_dict:
      return
    peaks = sum([len(ms2.peaks) if ms2.peaks is not None else 0 for ms2 in ms2_list])
    for (dictionary, tax) in [(self.kingdom_dict, 'kingdom'), (self.super_class_dict, 'super_class'), (self.sub_class_dict, 'sub_class'), (self.class_dict, 'class')]:
      tax_instance = taxonomy_dict.get(tax)
      if tax_instance:
        if dictionary.get(tax_instance):
          total_ms2, avg, count = dictionary.get(tax_instance)
          total_ms2 += 1
          if float(count + len(ms2_list)) == 0:
            avg = 0.
          else:
            avg = (avg*count + peaks)/float(count + len(ms2_list))
          count += len(ms2_list)
          dictionary[tax_instance] = [total_ms2, avg, count]
        else:
          if len(ms2_list):
            initial_avg = peaks/float(len(ms2_list))
          else:
            initial_avg = 0.
          dictionary[tax_instance] = [1., initial_avg, len(ms2_list)]



  def runStats(self):
    for metabolite in self.metabolite_dict.itervalues():
      self.updateIonization(metabolite.taxonomy_dict, metabolite.MS2)
      self.ms2_per_molecule.append(len(metabolite.MS2))
      for ms2 in metabolite.MS2:
        if ms2.peaks:
          self.peaks_per_ms2.append(len(ms2.peaks))
        else:
          self.peaks_per_ms2.append(0)
      self.updateCountingDicts(metabolite.taxonomy_dict, metabolite.MS2)


    # #graphing MS2/molecule
    # n, bins, patches = plt.hist(self.ms2_per_molecule, 40, facecolor='g', alpha=0.75)
    # plt.xlabel("MS2s/molecule")
    # plt.ylabel("Frequency of Occurrence")
    # plt.title('Exploring MS2s/molecule')
    # plt.axis([0, 40, 0, 500])
    # plt.grid(True)
    # plt.show()

    ##graphing peaks/MS2
    # n, bins, patches = plt.hist(self.peaks_per_ms2, 100, facecolor='g', alpha=0.75)
    # plt.xlabel("Peaks/MS2")
    # plt.ylabel("Frequency of Occurrence")
    # plt.title('Exploring Peaks/MS2')
    # plt.axis([0, 1500, 0, 100])  #[400, 1200, 4000]
    # plt.grid(True)
    # plt.annotate('heights: 4000, 1200, 400, 140', xy=(150, 75), xytext=(250, 75),
    #         arrowprops=dict(facecolor='black', shrink=0.05),
    #         )
    # plt.show()



    

 



    # ##graphing kingdom ionizations
    # dictionary = self.kingdom_ionization_dict
    # tax = 'kingdom'
    # values_array = np.zeros((1,3))
    # columns = ['tax', 'ionization', 'percentage']
    # df = pd.DataFrame(values_array, columns=columns)

    # for tax_instance, values in dictionary.iteritems():
    #   total = float(sum(values))
    #   for index in xrange(4):
    #     if index == 0:
    #       fixed = 0
    #     else:
    #       fixed = -1.*index
    #     df2 = pd.DataFrame(np.array([tax_instance, 1.*index, float(values[index])/total*100]).reshape(1,3), columns=columns)
    #     df = df.append(df2)
    # df = df.iloc[1:]
    # df['x']=sorted(range(df.shape[0]/4)*4)


    # print df
    # csv_path = os.getcwd() + '/ionization.csv'
    # df = pd.DataFrame.from_csv(csv_path).sort(columns='ionization')
    # print df
    # x = df.loc[:, 'x'].tolist()
    # tax = df.loc[:, 'tax'].tolist()[:len(x)/4]
    # ionization = df.loc[:, 'ionization'].tolist()
    # percentage = df.loc[:, 'percentage'].tolist()
    # tax_instances = len(x)/4


    # width = 0.35
    # colors = ['#d62728', '#0000ff', '#800080', '#00ff00']
    # graphs = []
    # for bar_graph_index in xrange(4):
    #   bar_x = np.arange(tax_instances)
    #   bar_height = tuple(percentage[tax_instances*bar_graph_index:tax_instances*(bar_graph_index+1)])
    #   print bar_height
    #   graph_color = colors[bar_graph_index]
    #   plot = plt.bar(bar_x, bar_height, width, color=graph_color)
    #   graphs.append(plot)
    # plt.ylabel('Percent')
    # plt.title('Ionization Type Percentages')
    # plt.xticks(bar_x, (tax))
    # plt.yticks(np.arange(0, 110, 10))
    # plt.ylim(0, 150)
    # plt.legend(tuple([graph[0] for graph in graphs]), ('Positive', 'Negative', 'N/A', 'None'))
    # print df
    # plt.show()



    # ###graphing super class ionizations
    #     dictionary = self.kingdom_ionization_dict
    # tax = 'kingdom'
    # values_array = np.zeros((1,3))
    # columns = ['tax', 'ionization', 'percentage']
    # df = pd.DataFrame(values_array, columns=columns)

    # for tax_instance, values in dictionary.iteritems():
    #   total = float(sum(values))
    #   for index in xrange(4):
    #     if index == 0:
    #       fixed = 0
    #     else:
    #       fixed = -1.*index
    #     df2 = pd.DataFrame(np.array([tax_instance, 1.*index, float(values[index])/total*100]).reshape(1,3), columns=columns)
    #     df = df.append(df2)
    # df = df.iloc[1:]
    # df['x']=sorted(range(df.shape[0]/4)*4)


    # # print df
    # csv_path = os.getcwd() + '/ionization.csv'
    # df = pd.DataFrame.from_csv(csv_path).sort(columns='ionization')
    # print df
    # x = df.loc[:, 'x'].tolist()
    # tax = df.loc[:, 'tax'].tolist()[:len(x)/4]
    # ionization = df.loc[:, 'ionization'].tolist()
    # percentage = df.loc[:, 'percentage'].tolist()
    # tax_instances = len(x)/4


    # width = 0.35
    # colors = ['#d62728', '#0000ff', '#800080', '#00ff00']
    # graphs = []
    # for bar_graph_index in xrange(4):
    #   bar_x = np.arange(tax_instances)
    #   bar_height = tuple(percentage[tax_instances*bar_graph_index:tax_instances*(bar_graph_index+1)])
    #   print bar_height
    #   graph_color = colors[bar_graph_index]
    #   plot = plt.bar(bar_x, bar_height, width, color=graph_color)
    #   graphs.append(plot)
    # plt.ylabel('Percent')
    # plt.title('Ionization Type Percentages')
    # plt.xticks(bar_x, (tax))
    # plt.yticks(np.arange(0, 110, 10))
    # plt.ylim(0, 150)
    # plt.legend(tuple([graph[0] for graph in graphs]), ('Positive', 'Negative', 'N/A', 'None'))
    # print df
    # plt.show()

    # assert False






    # x = df.loc[:, 'x'].tolist()
    # tax = df.loc[:, 'tax'].tolist()[:len(x)/4]
    # ionization = df.loc[:, 'ionization'].tolist()
    # percentage = df.loc[:, 'percentage'].tolist()
    # tax_instances = len(x)/4

    # width = 0.35
    # colors = ['#d62728', '#0000ff', '#800080', '#00ff00']
    # graphs = []
    # for bar_graph_index in xrange(4):
    #   bar_x = np.arange(tax_instances)
    #   bar_height = tuple(percentage[tax_instances*bar_graph_index:tax_instances*(bar_graph_index+1)])
    #   print bar_height
    #   graph_color = colors[bar_graph_index]
    #   plot = plt.bar(bar_x, bar_height, width, color=graph_color)
    #   graphs.append(plot)
    # plt.ylabel('Percent')
    # plt.title('Ionization Type Percentages')
    # plt.xticks(bar_x, (tax))
    # plt.yticks(np.arange(0, 110, 10))
    # plt.ylim(0, 150)
    # plt.legend(tuple([graph[0] for graph in graphs]), ('Positive', 'Negative', 'N/A', 'None'))
    # plt.show()






    # dictionary = self.super_class_ionization_dict
    # tax = 'super_class'
    # values_array = np.zeros((1,3))
    # columns = ['tax', 'ionization', 'percentage']
    # df = pd.DataFrame(values_array, columns=columns)

    # for tax_instance, values in dictionary.iteritems():
    #   total = float(sum(values))
    #   for index in xrange(4):
    #     if index == 0:
    #       fixed = 0
    #     else:
    #       fixed = -1.*index
    #     df2 = pd.DataFrame(np.array([tax_instance, 1.*index, float(values[index])/total*100]).reshape(1,3), columns=columns)
    #     df = df.append(df2)
    # df = df.iloc[1:]
    # df['x']=sorted(range(df.shape[0]/4)*4)



    # x = df.loc[:, 'x'].tolist()
    # tax = df.loc[:, 'tax'].tolist()[:len(x)/4]
    # ionization = df.loc[:, 'ionization'].tolist()
    # percentage = df.loc[:, 'percentage'].tolist()
    # tax_instances = len(x)/4

    # width = 0.35
    # colors = ['#d62728', '#0000ff', '#800080', '#00ff00']
    # graphs = []
    # for bar_graph_index in xrange(4):
    #   bar_x = np.arange(tax_instances)
    #   bar_height = tuple(percentage[tax_instances*bar_graph_index:tax_instances*(bar_graph_index+1)])
    #   print bar_height
    #   graph_color = colors[bar_graph_index]
    #   plot = plt.bar(bar_x, bar_height, width, color=graph_color)
    #   graphs.append(plot)
    # plt.ylabel('Percent')
    # plt.title('Ionization Type Percentages')
    # plt.xticks(bar_x, (tax))
    # plt.yticks(np.arange(0, 110, 10))
    # plt.ylim(0, 150)
    # plt.legend(tuple([graph[0] for graph in graphs]), ('Positive', 'Negative', 'N/A', 'None'))
    # plt.show()



    # assert False


    for (tax, dictionary) in [('kingdom', self.kingdom_ionization_dict), ('super_class', self.super_class_ionization_dict), ('class', self.class_ionization_dict), ('sub_class', self.sub_class_ionization_dict)]:
      values_array = np.zeros((1,3))
      columns = ['tax', 'ionization', 'percentage']
      df = pd.DataFrame(values_array, columns=columns)
      for tax_instance, values in dictionary.iteritems():
        total = float(sum(values))
        for index in xrange(4):
          if index == 0:
            fixed = 0
          else:
            fixed = -1.*index
          df2 = pd.DataFrame(np.array([tax_instance, 1.*index, float(values[index])/total*100]).reshape(1,3), columns=columns)
          df = df.append(df2)
      df = df.iloc[1:]
      df['x']=sorted(range(df.shape[0]/4)*4)
      df.to_csv(os.getcwd() + '/ionization.csv')


      csv_path = os.getcwd() + '/ionization.csv'
      df = pd.DataFrame.from_csv(csv_path).sort(columns='ionization')
      x = df.loc[:, 'x'].tolist()
      tax = df.loc[:, 'tax'].tolist()[:len(x)/4]
      ionization = df.loc[:, 'ionization'].tolist()
      percentage = df.loc[:, 'percentage'].tolist()
      tax_instances = len(x)/4


      width = 0.35
      colors = ['#d62728', '#0000ff', '#800080', '#00ff00']
      graphs = []
      for bar_graph_index in xrange(4):
        bar_x = np.arange(tax_instances)
        bar_height = tuple(percentage[tax_instances*bar_graph_index:tax_instances*(bar_graph_index+1)])
        print bar_height
        graph_color = colors[bar_graph_index]
        plot = plt.bar(bar_x, bar_height, width, color=graph_color)
        graphs.append(plot)
      plt.ylabel('Percent')
      plt.title('Ionization Type Percentages')
      plt.xticks(bar_x, (tax))
      plt.yticks(np.arange(0, 110, 10))
      plt.ylim(0, 150)
      plt.legend(tuple([graph[0] for graph in graphs]), ('Positive', 'Negative', 'N/A', 'None'))
      plt.xticks(rotation=300)
      print df
      plt.show()














      # sns.set(style="whitegrid")
      ####CLAIRE FIX THIS LINE PLEASE######
      # sns.factorplot(kind='bar', x='x', y='percentage', hue_order=['3', '2', '1', '0'], hue='ionization', data=df)
      ## sns.factorplot(kind='bar', x='x', y='percentage', hue_order=['Positive', 'Negative', 'N/A', 'None'], hue='ionization', data=df)
      # plt.show()


    #graphing total MS2s and avg peaks/MS2 at each taxonomy level
    #one column with the specific taxonomy_instance
    #a column with the x values
    #a column with the total MS2's for that particular taxonomy instance
    #a column with the average peaks/MS2 at that particular taxonomy instance
    #a column that provides hues
    assert False
    for (tax, dictionary) in [('kingdom', self.kingdom_dict), ('super_class', self.super_class_dict), ('class', self.class_dict), ('sub_class', self.sub_class_dict)]:
      values_array = np.zeros((1,4))
      columns = ['tax', 'MS2s', 'avg peaks/MS2', 'x']
      df = pd.DataFrame(values_array, columns=columns)
      counter = 0
      for tax_instance, (MS2s, peaks_per, _) in dictionary.iteritems():
        df2 = pd.DataFrame(np.array([tax_instance, MS2s, peaks_per, counter]).reshape(1,4), columns=columns)
        df = df.append(df2)
        counter += 1
      df = df.iloc[1:]
      sns.set(style="whitegrid")
      ###CLAIRE FIX THIS PLEASE####
      sns.factorplot(kind='bar', x='x', y='MS2s', hue='x', data=df)
      plt.title('Total number of MS2s for each %s' % (tax))
      plt.ylabel('Number of molecules')
      plt.show()

      ###CLAIRE FIX THIS PLEASE####
      sns.factorplot(kind='bar', x='x', y='avg peaks/MS2', data=df)
      plt.title('Avg number of peaks per MS2 for each %s' % (tax))
      plt.ylabel('Number of molecules')
      plt.show()



# def countMS2AndMetabolites():
#   csv_path = os.getcwd() + '/../../data/clean/metabolites_and_spectra.csv'
#   csv_file = open(csv_path, 'r')
#   reader = csv.reader(csv_file)
#   metabolite_count = 0
#   ms2_count = 0
#   for row in reader:
#     if row[0]:
#       metabolite_count += 1
#     else:
#       ms2_count += 1

#   print metabolite_count, 'metabolites found'
#   print ms2_count, 'MS2s found'




if __name__=='__main__':
  stats = StatsMaker()
  stats.runStats()

  # csv_path = os.getcwd() + '/ionization.csv'
  # df = pd.DataFrame.from_csv(csv_path).sort(columns='ionization')
  # x = df.loc[:, 'x'].tolist()
  # tax = df.loc[:, 'tax'].tolist()[:len(x)/4]
  # ionization = df.loc[:, 'ionization'].tolist()
  # percentage = df.loc[:, 'percentage'].tolist()
  # tax_instances = len(x)/4


  # width = 0.35
  # colors = ['#d62728', '#0000ff', '#800080', '#00ff00']
  # graphs = []
  # for bar_graph_index in xrange(4):
  #   bar_x = np.arange(tax_instances)
  #   bar_height = tuple(percentage[tax_instances*bar_graph_index:tax_instances*(bar_graph_index+1)])
  #   print bar_height
  #   graph_color = colors[bar_graph_index]
  #   plot = plt.bar(bar_x, bar_height, width, color=graph_color)
  #   graphs.append(plot)
  # plt.ylabel('Percent')
  # plt.title('Ionization Type Percentages')
  # plt.xticks(bar_x, (tax))
  # plt.yticks(np.arange(0, 110, 10))
  # plt.ylim(0, 150)
  # plt.legend(tuple([graph[0] for graph in graphs]), ('Positive', 'Negative', 'N/A', 'None'))
  # print df
  # plt.show()








#  N = 2
# menMeans = (20, 35)
# womenMeans = (25, 32)
# ind = np.arange(N)    # the x locations for the groups
# width = 0.35       # the width of the bars: can also be len(x) sequence

# p1 = plt.bar(ind, menMeans, width, color='#d62728')
# p2 = plt.bar(ind, womenMeans, width,
#              bottom=menMeans)

# plt.ylabel('Scores')
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# # plt.legend((p1[0], p2[0]), ('Men', 'Women'))

# plt.show()







