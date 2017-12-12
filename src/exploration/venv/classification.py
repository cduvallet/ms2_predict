import sys
import numpy as np
import pdb

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
import pandas as pd
import os

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split

import datetime
import csv


class Trainer:
  def __init__(self):
    self.classes = ['Organic acids and derivatives;', 'Lipids and lipid-like molecules;', 'Organoheterocyclic compounds;', 'Organic oxygen compounds;', 'Benzenoids;', 'Nucleosides, nucleotides, and analogues;', 'Carboxylic acids and derivatives;', 'Steroids and steroid derivatives;', 'Organooxygen compounds;', 'Fatty Acyls;', 'Benzene and substituted derivatives;', 'Phenylpropanoids and polyketides;']
    self.subclasses = ['carboxylic acids and derivatives', 'fatty acyls', 'organooxygen compounds', 'amino acids, peptides, and analogues', 'benzene and substituted derivatives', 'imidazopyrimidines', 'fatty acids and conjugates', 'phenols', 'indoles and derivatives', 'carbohydrates and carbohydrate conjugates', 'organonitrogen compounds', 'Fatty acids and conjugates', 'Phenols', 'Indoles and derivatives', 'Carbohydrates and carbohydrate conjugates', 'Organonitrogen compounds']
    main_csv = open(os.getcwd() + '/main.csv', 'w+')
    self.writer = csv.writer(main_csv)
    self.id_counter = 1
    self.writer.writerow(['ID', 'Type', 'C', 'gamma', 'n_estimators', 'max_features', 'min_samples_leaf', 'val accuracy', 'test accuracy', 'feature importances'])



  def preprocess(self, path, tax_type):
    #last four rows are inchi, kingdom, sub_class, class
    df = pd.read_table(path)
    df = df.rename(index=str, columns={"class": "_class"})


    filter_df = df
    if tax_type == 'subclass':
      pdb.set_trace()
      filter_df = filter_df[filter_df['sub_class'].isin(self.subclasses)]
      if 'ms2lda' in path:
        self.X = filter_df.loc[:, 'motif_0-overlap':'motif_99-prob']
      else:
        self.X = filter_df.loc[:, '0':'1980']
      self.Y = filter_df.loc[:, 'sub_class']

    elif tax_type == 'class':
      filter_df = filter_df[filter_df['_class'].isin(self.classes)]
      if 'ms2lda' in path:
        self.X = filter_df.loc[:, 'motif_0-overlap':'motif_99-prob']
      else:
        self.X = filter_df.loc[:, '0':'1980']
      self.Y = filter_df.loc[:, '_class']

    else:
      if 'ms2lda' in path:
        self.X = df.loc[:, 'motif_0-overlap':'motif_99-prob']
      else:
        self.X = df.loc[:, '0':'1980']
      self.Y = df.loc[:, 'kingdom']


  def partitionData(self):
    cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2)
    split = [ _ for _ in cv.split(self.X, self.Y)]

    self.X_val = pd.DataFrame(self.X.iloc[[split[0][0][0]]])
    self.Y_val = list(self.Y.iloc[[split[0][0][0]]])
    for index in split[0][0][1:]:
      val_row = self.X.iloc[[index]]
      self.X_val = self.X_val.append(val_row)
      val_label = self.Y.iloc[[index]]
      self.Y_val.append(val_label[0])
    self.Y_val = pd.DataFrame(self.Y_val)

    self.X_test = pd.DataFrame(self.X.iloc[[split[0][1][0]]])
    self.Y_test = list(self.Y.iloc[[split[0][1][0]]])
    for index in split[0][0][1:]:
      test_row = self.X.iloc[[index]]
      self.X_test = self.X_test.append(test_row)
      test_label = self.Y.iloc[[index]]
      self.Y_test.append(test_label[0])
    self.Y_test = pd.DataFrame(self.Y_test)


  '''
  ID#.csv is the cv_results_ dictionary from k-folds validation
  main.csv records optimal parameters
  '''
  def tuning(self, tax_type):
    print 'Tuning SVM'
    C_range = np.logspace(-3, 3, 10)
    gamma_range = np.logspace(-9, 3, 10)
    # C_range = [1]
    # gamma_range = [1]
    param_grid = dict(gamma=gamma_range, C=C_range)
    cv = StratifiedShuffleSplit(n_splits=4, test_size=0.2)
    grid = GridSearchCV(SVC(kernel='rbf', tol=1e-3, decision_function_shape='ovr'), param_grid=param_grid, cv=cv, return_train_score=True)

    grid.fit(self.X_val, self.Y_val)
    c_opt = grid.best_params_['C']
    gamma_opt = grid.best_params_['gamma']
    # score_opt = grid.best_score_
    pd.DataFrame(grid.cv_results_).to_csv(str(self.id_counter)+'.csv')

    #retrain on validation data
    clf = SVC(kernel='rbf', tol=1e-3, decision_function_shape='ovr', C=c_opt, gamma=gamma_opt, probability=True)
    clf.fit(self.X_val, self.Y_val)

    val_predictions = clf.predict(self.X_val)
    val_probs = clf.predict_proba(self.X_val)
    val_score = clf.score(self.X_val, self.Y_val)


    Y_val_storage = self.Y_val
    Y_val_storage = Y_val_storage.rename(index=str, columns={0: 'Label'})
    Y_val_storage['Predicted'] = val_predictions

    for _class in xrange(val_probs.shape[1]):
      Y_val_storage['prob' + str(_class+1)] = val_probs[:, _class]
    file_name = str(self.id_counter) + '_val_data.csv'
    Y_val_storage.to_csv(file_name)


    test_predictions = clf.predict(self.X_test)
    test_probs = clf.predict_proba(self.X_test)
    test_score = clf.score(self.X_test, self.Y_test)

    Y_test_storage = self.Y_test
    Y_test_storage = Y_test_storage.rename(index=str, columns={0: 'Label'})
    Y_test_storage['Predicted'] = test_predictions

    for _class in xrange(test_probs.shape[1]):
      Y_test_storage['prob' + str(_class+1)] = test_probs[:, _class]
    file_name = str(self.id_counter) + '_test_data.csv'
    Y_test_storage.to_csv(file_name)


    #confusion matrices
    test_cm = confusion_matrix(self.Y_test, test_predictions)
    np.savetxt(str(self.id_counter) + '_test_cm.txt', test_cm)
    val_cm = confusion_matrix(self.Y_val, val_predictions)
    np.savetxt(str(self.id_counter) + '_val_cm.txt', test_cm)

    # self.AUC(tax_type)
    
    row = [self.id_counter, 'SVM', c_opt, gamma_opt, '', '', '', val_score, test_score]
    self.writer.writerow(row)

    self.id_counter += 1


    print 'Tuning Random Forest'
    max_features = ['auto', 'log2']
    # max_features = ['auto']
    min_samples_leaf = [1, 2, 3]
    # min_samples_leaf = [1]
    n_estimators = [100, 1000, 5000]
    # n_estimators = [10]
    param_grid = dict(max_features=max_features, min_samples_leaf=min_samples_leaf, n_estimators=n_estimators)
    cv = StratifiedShuffleSplit(n_splits=4, test_size=0.2)
    grid = GridSearchCV(RandomForestClassifier(), param_grid=param_grid, cv=cv)
    grid.fit(self.X, self.Y)
    max_features_opt = grid.best_params_['max_features']
    min_samples_leaf_opt = grid.best_params_['min_samples_leaf']
    n_estimators_opt = grid.best_params_['n_estimators']
    pd.DataFrame(grid.cv_results_).to_csv(str(self.id_counter)+'.csv')

    #retrain on validation data
    clf = RandomForestClassifier(max_features=max_features_opt, min_samples_leaf=min_samples_leaf_opt, n_estimators=n_estimators_opt)
    clf.fit(self.X_val, self.Y_val)

    val_predictions = clf.predict(self.X_val)
    val_probs = clf.predict_proba(self.X_val)
    val_score = clf.score(self.X_val, self.Y_val)

    Y_val_storage = self.Y_val
    Y_val_storage = Y_val_storage.rename(index=str, columns={0: 'Label'})
    Y_val_storage['Predicted'] = val_predictions

    for _class in xrange(val_probs.shape[1]):
      Y_val_storage['prob' + str(_class+1)] = val_probs[:, _class]
    file_name = str(self.id_counter) + '_val_data.csv'
    Y_val_storage.to_csv(file_name)

    test_predictions = clf.predict(self.X_test)
    test_probs = clf.predict_proba(self.X_test)
    test_score = clf.score(self.X_test, self.Y_test)

    Y_test_storage = self.Y_test
    Y_test_storage = Y_test_storage.rename(index=str, columns={0: 'Label'})
    Y_test_storage['Predicted'] = test_predictions

    for _class in xrange(test_probs.shape[1]):
      Y_test_storage['prob' + str(_class+1)] = test_probs[:, _class]
    file_name = str(self.id_counter) + '_test_data.csv'
    Y_test_storage.to_csv(file_name)

    #confusion matrices
    test_cm = confusion_matrix(self.Y_test, test_predictions)
    np.savetxt(str(self.id_counter) + '_test_cm.txt', test_cm)
    val_cm = confusion_matrix(self.Y_val, val_predictions)
    np.savetxt(str(self.id_counter) + '_val_cm.txt', test_cm)

    # self.AUC(tax_type)
    
    row = [self.id_counter, 'RF', '', '', n_estimators_opt, max_features_opt, min_samples_leaf_opt, val_score, test_score, str(clf.feature_importances_)]
    self.writer.writerow(row)

    self.id_counter += 1


  def AUC(self, tax_type):
    #calculate AUC
    #must transition to 'binary' classification
    if tax_type == 'kingdom': 
      classes = 2
    elif tax_type == 'subclass': classes = len(self.subclasses)
    elif tax_type == 'class': classes = len(self.classes)

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    distances = clf.decision_function(self.X_test)
    for i in range(classes):
      fpr[i], tpr[i], _ = roc_curve(np.array(self.Y_test.as_matrix())[:, i], distances.reshape((distances.shape[0], 1))[:, i], pos_label='Chemical entities')
      roc_auc[i] = auc(fpr[i], tpr[i])

    #plot AUC curve
    plt.figure()
    lw = 2
    plt.plot(fpr[2], tpr[2], color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()




  def main(self):
    binned_path = os.getcwd() + '/../../../data/feature_tables/mz_integer.all_scans.txt'
    ms2lda_path = os.getcwd() + '/../../../data/feature_tables/ms2lda_feature_table.merged_spectra.txt'


    for path in (binned_path, ms2lda_path):
    # for path in [binned_path]:
      for tax_type in ('kingdom', 'subclass', 'class'):
      # for tax_type in ['subclass']:
        print 'Running pipeline for data=%s, tax_type=%s' % (['binned data', 'ms2lda data'][[binned_path, ms2lda_path].index(path)], tax_type)
        self.preprocess(path, tax_type)
        self.partitionData()
        self.tuning(tax_type)



if __name__ == '__main__':
  trainer = Trainer()
  trainer.main()


'''
Split the data 80/20
4-fold cross validation to get optimal parameters
store parameters, cv thing
retrain on all validation data, store feature importances
apply classifier to test data. Store label, predicted label, probabilities
save percent accuracy, confusion matrix, AUC for validation and training data
'''


