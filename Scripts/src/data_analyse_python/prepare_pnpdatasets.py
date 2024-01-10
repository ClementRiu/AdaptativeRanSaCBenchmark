import numpy as np
import os
import random
from scipy.spatial.transform import Rotation as R

class camera:
    def __init__(self, data_):
        self._cam_id = int(data_[0])
        self._cam_model = data_[1]

        self._cam_dim = np.array([data_[2], data_[3]]).astype(float)

        if self._cam_model == "SIMPLE_RADIAL":
            self._cam_param = np.array([data_[4], 0, data_[5], data_[4], data_[6]]).astype(float)
            self._distortion_param = float(data_[7])
        else:
            print("NO MODEL AVAILABLE.")

class image:
    def __init__(self, data_, points2d_, image_path):
        self._im_id = int(data_[0])

        self._im_quaternion = np.array([data_[1], data_[2], data_[3], data_[4]]).astype(float)
        self._im_translation = np.array([data_[5], data_[6], data_[7]]).astype(float)

        self._cam_id = int(data_[8])
        self._im_name = image_path + data_[9]

        self._points2d = {-1: []}
        for index in range(0, len(points2d_), 3):
            pt3d_id = int(points2d_[index + 2])
            pt2d = np.array([points2d_[index], points2d_[index + 1]]).astype(float)
            if pt3d_id == -1:
                self._points2d[-1].append(pt2d)
            else:
                self._points2d[pt3d_id] = pt2d

    def getMat(self):
        self._R = np.zeros((3, 3))
        qr, qi, qj, qk = self._im_quaternion
        s = 1 / self._im_quaternion.dot(self._im_quaternion)
        self._R[0, 0] = 1 - 2 * s * (qj * qj + qk * qk)
        self._R[0, 1] = 2 * s * (qi * qj - qk * qr)
        self._R[0, 2] = 2 * s * (qi * qk + qj * qr)
        self._R[1, 0] = 2 * s * (qi * qj + qk * qr)
        self._R[1, 1] = 1 - 2 * s * (qi * qi + qk * qk)
        self._R[1, 2] = 2 * s * (qj * qk - qi * qr)
        self._R[2, 0] = 2 * s * (qi * qk - qj * qr)
        self._R[2, 1] = 2 * s * (qj * qk + qi * qr)
        self._R[2, 2] = 1 - 2 * s * (qi * qi + qj * qj)

    def addCamera(self, camera):
        self._camera = camera

    def correct_radial(self):
        self._points2d_corr = {-1: []}
        for key, point2D in self._points2d.items():
            if key != -1 :
                self._points2d_corr[key] = self.apply_correction(point2D)
                # self._points2d_corr[key] = point2D

    def apply_correction(self, point2D):
        point_ref = np.array([self._camera._cam_param[2], self._camera._cam_param[4]])
        focal = self._camera._cam_param[0]
        point2D -= point_ref
        point2D /= focal
        r = (point2D).dot(point2D)
        point2D += point2D * self._camera._distortion_param * r
        point2D = point2D * focal + point_ref
        if point2D[0] >= 0 and point2D[0] < self._camera._cam_dim[0] and point2D[1] >= 0 and point2D[1] < self._camera._cam_dim[1]:
            return point2D

    def add3Dpts(self, pts3d_):
        for pt3D_id in self._points2d_corr.keys():
            if pt3D_id != -1:
                    self._points2d_corr[pt3D_id] = np.concatenate((self._points2d_corr[pt3D_id], pts3d_[pt3D_id]._pt3d_loc))

    def project3D2D(self):
        for pt3D_id in self._points2d_corr.keys():
            if pt3D_id != -1:
                point3D = self._points2d_corr[pt3D_id][2:]
                xyz = self._R.dot(point3D) + self._im_translation.T
                u = xyz[0] * self._camera._cam_param[0] + xyz[1] * self._camera._cam_param[1] + xyz[2] * self._camera._cam_param[2]
                v = xyz[1] * self._camera._cam_param[3] + xyz[2] * self._camera._cam_param[4]
                if xyz[2] > 0:
                    u /= xyz[2]
                    v /= xyz[2]
                    if u >= 0 and u < self._camera._cam_dim[0] and v >= 0 and v < self._camera._cam_dim[1]:
                    # print(self._points2d_corr[pt3D_id][:2], u, v)
                        self._points2d_corr[pt3D_id][0] = u
                        self._points2d_corr[pt3D_id][1] = v

    def saveData(self, path_calib, path_points, path_image, path_RT):
        np.savetxt(path_calib, self._camera._cam_param)

        allPoints = np.array([
            pt2D3D for pt3D_id, pt2D3D in self._points2d_corr.items() if pt3D_id != -1
        ])
        np.savetxt(path_points, allPoints)

        os.system("cp {} {}".format(self._im_name, path_image))

        RT_values = np.concatenate((self._R, [self._im_translation]))

        np.savetxt(path_RT, RT_values)

class point:
    def __init__(self, data_):
        self._pt3d_id = int(data_[0])

        self._pt3d_loc = np.array([data_[1], data_[2], data_[3]]).astype(float)
        self._pt3d_rgb = np.array([data_[4], data_[5], data_[6]]).astype(float)
        self._pt3d_err = float(data_[7])

        self._track2d = []
        for index in range(8, len(data_), 2):
            self._track2d.append((int(data_[index]), int(data_[index + 1])))


if __name__=="__main__":
    dataset_path = "/home/riuclement//Documents/datasets/megadepth_sfm/0009/"
    image_path = dataset_path + "images/"
    data_path = "sparse/manhattan/0/"
    cameras_path = dataset_path + data_path + "cameras.txt"
    images_path = dataset_path + data_path + "images.txt"
    points_path = dataset_path + data_path + "points3D.txt"

    output_path = "/home/riuclement//Documents/datasets/megadepth_sfm/pnp_ready/{}_{}/".format("0009", "{}")
    output_names = ["calib_matrix.txt", "2D3D_pts.txt", "im.jpg", "RT.txt"]

    cameras = []
    with open(cameras_path, 'r') as cameras_file:
        cameras = cameras_file.readlines()

    images = []
    with open(images_path, 'r') as images_file:
        images = images_file.readlines()

    points = []
    with open(points_path, 'r') as points_file:
        points = points_file.readlines()

    cameras = [elem.strip().split(' ') for elem in cameras]
    images = [elem.strip().split(' ') for elem in images]
    points = [elem.strip().split(' ') for elem in points]

    cameras = [camera(elem) for elem in cameras[3:]]
    images = [image(images[index], images[index + 1], image_path) for index in range(4, len(images), 2)]
    points = [point(elem) for elem in points[4:]]

    cameras = {elem._cam_id : elem for elem in cameras}
    points = {elem._pt3d_id : elem for elem in points}

    selected_images = random.choices(images, k=2)

    for image in selected_images:
        image.addCamera(cameras[image._cam_id])
        image.getMat()
        image.correct_radial()
        image.add3Dpts(points)
        # image.project3D2D()
        os.system("mkdir -p {}".format(output_path.format(image._im_id)))
        image.saveData(
            output_path.format(image._im_id) + output_names[0],
            output_path.format(image._im_id) + output_names[1],
            output_path.format(image._im_id) + output_names[2],
            output_path.format(image._im_id) + output_names[3],
        )
