#!/usr/bin/env python
# coding: utf-8

from math import inf
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score


def build_AM(x,y):
    """Function to build the A matrix.

    Keyword arguments:
    x -- the input data numpy array in the form nxm (samplesxfeatures).
    y -- the numpy array that represents the classes for each sample.
    """
    Am = np.zeros((x.shape[1],len(np.unique(y))))
    for feat in range(x.shape[1]):
        for lab in np.unique(y):
            x_fit = x[np.where(y == lab)[0],feat].reshape(-1, 1)
            params = {'bandwidth': np.linspace(0.01, 3, 30)}
            grid = GridSearchCV(KernelDensity(), params, cv=5)
            grid.fit(x_fit)
            kde = grid.best_estimator_
            #h = np.std(x_fit)*(4/3/len(x_fit))**(1/5)
            #print(h)
            #kde = KernelDensity(bandwidth=max(h,0.5)).fit(x_fit)
            Am[feat,lab] = x_fit[np.argmax(np.exp(kde.score_samples(x_fit)))]
    return Am

def new_sample_prediction(Am, b_mat, pred=None):
    """Solves the Ax=b linear system and returns the class predicted.

    Keyword arguments:
    Am -- The A matrix.
    b_mat -- The b matrix representing the sample(s) from which we want predictions to.
    pred -- If None the entire vector is returned. Otherwise may be one of the following strings: "HV" so the prediction returned is the Highest Value, 
    "HM" for the Highest Magnitude, or "LM_1" for the Lowest Magnitude closest to 1.
    """
    assert Am.shape[0] == b_mat.shape[0], "The number of features on b (first dimension) doesn't match the first dimension from A."
    output = list(np.linalg.lstsq(Am, b_mat,rcond=None))
    output[0] = np.transpose(output[0])
    if not pred:
        return output
    elif pred == "HV":
        return np.argmax(output[0], axis=1)
    elif pred == "HM":
        return np.argmax(np.absolute(output[0]), axis=1)
    elif pred == "LM_1":
        return np.argmin(np.absolute(output[0]-1), axis=1)

def new_sample_prediction_lu(Am, b_mat, pred=None):
    """Solves the Ax=b linear system using LU factorization returns the class predicted.

    Keyword arguments:
    Am -- The A matrix.
    b_mat -- The b matrix representing the sample(s) from which we want predictions to.
    pred -- If None the entire vector is returned. Otherwise may be one of the following strings: "HV" so the prediction returned is the Highest Value, 
    "HM" for the Highest Magnitude, or "LM_1" for the Lowest Magnitude closest to 1.
    """
    assert Am.shape[0] == b_mat.shape[0], "The number of features on b (first dimension) doesn't match the first dimension from A."
    lu, piv = lu_factor(Am)
    output = list(lu_solve((lu, piv), b_mat))
    output = np.transpose(output)
    if not pred:
        return output
    elif pred == "HV":
        return np.argmax(output, axis=1)
    elif pred == "HM":
        return np.argmax(np.absolute(output), axis=1)
    elif pred == "LM_1":
        return np.argmin(np.absolute(output+1), axis=1)


def confusion_matrix_print(real, pred, labels):
    cmsk = confusion_matrix(list(real), list(pred))
    base = np.unique(real)
    base.sort()
    print("{:12}".format(""),end="")
    for e in list(labels.inv.values()):
        print("{:^13}".format(str(e)+"P"),end="")
    for r in base:
        print("\n{:>12}".format(str(labels.inv[r])+"A"),end="")
        for p in base:
            print("{:^13}".format(str(cmsk[r][p])),end="")
    print("\n\nAccuracy: {:.2f}%\n{}/{} samples classified correctly".format(
        100*accuracy_score(list(real), list(pred)), 
        accuracy_score(list(real), list(pred), normalize=False),
        len(list(real))))
    return 100*accuracy_score(list(real), list(pred))


def fbyc_plot(x_train, y_train, x_test, y_test, y_test_pred, labels_dict, feature_names, plot_correct = False, plot_incorrect = False, save_dir="images", filename="CbyF_incorrect_samples"):
    from plotly import tools
    import plotly.plotly as py
    import plotly.graph_objs as go
    from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
    import plotly.io as pio
    import plotly.figure_factory as ff
    import colorlover as cl
    #Uncomment below if using jupyter notebook
    #init_notebook_mode(connected=True)
    #compute x_lim
    x_lim_feature = dict()
    for feat in range(x_train.shape[1]):
        x_lim_feature[feat] = [inf, -inf]
        for key, label in labels_dict.items():
            if len(np.unique(x_train[np.where(y_train == label)[0],feat])) > 1:
                x_lim_feature[feat] =  [min(x_lim_feature[feat][0], min(np.unique(np.vstack((x_train, x_test))[np.where(np.hstack((y_train, y_test)) == label)[0],feat]))),
                                        max(x_lim_feature[feat][1], max(np.unique(np.vstack((x_train, x_test))[np.where(np.hstack((y_train, y_test)) == label)[0],feat])))]
    for key, value in x_lim_feature.items():
        diff = (value[1] - value[0])/10
        x_lim_feature[key][0] -= diff
        x_lim_feature[key][1] += diff
    #end compute x_lim
    title_list = ['F_{}'.format(i) for i in feature_names]
    title_list = title_list*len(labels_dict)
    fig = tools.make_subplots(rows=len(labels_dict), cols=x_train.shape[1], print_grid=False, subplot_titles=tuple(title_list))#, specs=[[{}]*30]*7, horizontal_spacing=1)
    k=1
    dict_label_plot = {p: True for p in labels_dict.keys()}
    plot_colors = cl.scales[str(len(labels_dict))]["div"]["RdYlBu"]
    for key_index, (key, label) in enumerate(labels_dict.items()):
        forceY = True
        for feat in range(x_train.shape[1]):
            hist_data = []
            group_labels = []
            hist_data.append(x_train[np.where(y_train == label)[0],feat].tolist())
            group_labels.append(label)
            if len(np.unique(x_train[np.where(y_train == label)[0],feat])) > 1:
                try:
                    trace1 = ff.create_distplot(hist_data, group_labels, show_hist=False, bin_size=0.05)["data"]
                    trace1[0]['showlegend'] = False
                    trace1[1]['showlegend'] = False
                    fig.append_trace(trace1[0], key_index+1, feat+1)
                    #Add Correct predictions to plot
                    correct_pre_ind = np.intersect1d(np.where(y_test == label)[0], np.where(y_test == y_test_pred)[0])
                    if len(correct_pre_ind) > 0 and plot_correct:
                        for correct_class in np.unique(y_test_pred[correct_pre_ind]):
                            correct_pre_ind = np.intersect1d(np.intersect1d(np.where(y_test == label)[0], 
                                                               np.where(y_test == y_test_pred)[0]),
                                                               np.where(correct_class == y_test_pred)[0])
                            x_train_f = x_test[correct_pre_ind,feat]
                            trace2 = go.Scatter(
                            x = x_train_f,
                            y = np.zeros_like(x_train_f),
                            mode = 'markers',
                            name = labels_dict.inv[correct_class],
                            marker={'color': plot_colors[correct_class], 'symbol': 4, 'size': 10},
                            showlegend = dict_label_plot[labels_dict.inv[correct_class]])
                            dict_label_plot[labels_dict.inv[correct_class]] = False
                            fig.append_trace(trace2, key_index+1, feat+1)
                    #Add Failed prediction to plot
                    incorrect_pre_ind = np.intersect1d(np.where(y_test == label)[0], np.where(y_test != y_test_pred)[0])
                    if len(incorrect_pre_ind) > 0 and plot_incorrect:
                        for incorrect_class in np.unique(y_test_pred[incorrect_pre_ind]):
                            incorrect_pre_ind = np.intersect1d(np.intersect1d(np.where(y_test == label)[0], 
                                                               np.where(y_test != y_test_pred)[0]),
                                                               np.where(incorrect_class == y_test_pred)[0])
                            x_train_f = x_test[incorrect_pre_ind,feat]
                            trace2 = go.Scatter(
                            x = x_train_f,
                            y = np.zeros_like(x_train_f),
                            mode = 'markers',
                            name = labels_dict.inv[incorrect_class],
                            marker={'color': plot_colors[incorrect_class], 'symbol': 4, 'size': 10},
                            showlegend = dict_label_plot[labels_dict.inv[incorrect_class]])
                            dict_label_plot[labels_dict.inv[incorrect_class]] = False
                            fig.append_trace(trace2, key_index+1, feat+1)
                    #Finish adding layout details
                    fig['layout']["xaxis{}".format(k)].update(zeroline = False, range = [x_lim_feature[feat][0],x_lim_feature[feat][1]])
                    if forceY:
                        fig['layout']['yaxis{}'.format(k)].update(title='{}'.format(key))
                        forceY = False
                except np.linalg.LinAlgError:
                    print("Feature {}, Class {} throws a LinAlgError. Let's not look at that...".format(feat+1, label))
            k+=1
    fig['layout'].update(height=1150, width=2250)
    iplot(fig)
    fig['layout'].update(height=2100, width=4500)
    pio.write_image(fig, '{}/{}.pdf'.format(save_dir, filename), scale=5)
