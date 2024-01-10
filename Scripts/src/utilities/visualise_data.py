import matplotlib.pyplot as plt
import numpy as np
import sys


def transform(model, x):
    return model.dot(np.array([x, 1]))

if __name__ == '__main__':

    # Loading of the file.
    with open(sys.argv[1]) as file:
        lines = file.readlines()
        char = lines[0]
        model = lines[1]
        inliers_index = lines[2]
        data = lines[3]

    width, height = char.strip().split(' ')
    width = int(width)
    height = int(height)

    # Loading of the model.
    model = np.array([float(value) for value in model.strip().split(' ') if value not in ['[', ']']])

    # Loading of the inlier index.
    inliers_index = np.array([int(value) for value in inliers_index.strip().split(' ')])

    # Loading of the two dimensional data.
    data1, data2 = data.split(';  ')
    data1 = [value for value in data1.split(' ') if value not in ['[', ']', '\n']]
    data2 = [value for value in data2.split(' ') if value not in ['[', ']', '\n']]
    data = np.array([[int(value1), int(value2)] for value1, value2 in zip(data1, data2)])

    # Loading of the inliers.
    inliers = data[inliers_index, :]
    # Loading of the outliers.
    if len(inliers_index) < 500:
        outliers = data[np.array([index for index in range(0, data.shape[0]) if index not in inliers_index]), :]
    else:
        outliers = np.array([[-1, -1]])

    print('The estimated model is: y = {} * x + {}.'.format(model[0], model[1]))
    print('There are {} points including {} inliers'.format(data.shape[0], inliers.shape[0]))

    # Creation of the estimated model line points.
    x0 = 0
    x1 = width
    y0 = transform(model, x0)
    y1 = transform(model, x1)

    plt.figure()
    plt.title('model y = {} * x + {}.'.format(model[0], model[1]))

    plt.xlim(0, width)
    plt.ylim(0, height)
    plt.xlabel('x')
    plt.ylabel('y')

    outplt = plt.scatter(outliers[:, 0], outliers[:, 1], s=30, c='r', alpha=0.5, marker='+', label='Outliers')
    inplt = plt.scatter(inliers[:, 0], inliers[:, 1], s=20, c='b', alpha=0.5, marker='+', label='Inliers')

    modplt, = plt.plot([x0, x1], [y0, y1], c='k', label='Estimated model')
    plt.legend(handles=[outplt, inplt, modplt], loc=(1.04, 0))
    plt.subplots_adjust(right=0.7)
    plt.show()
