import timeit
import skimage
import cv2


def belapsed_skimage(filename, *, number):
    t_read = timeit.timeit(lambda: skimage.io.imread(filename), number=number)
    img = skimage.io.imread(filename)
    t_write = timeit.timeit(
        lambda: skimage.io.imsave(filename, img), number=number)
    return 1000*t_read/number, 1000*t_write/number  # ms


def belapsed_cv2(filename, flag, *, number):
    t_read = timeit.timeit(lambda: cv2.imread(filename, flag), number=number)
    img = cv2.imread(filename, flag)
    t_write = timeit.timeit(lambda: cv2.imwrite(filename, img), number=number)
    return 1000*t_read/number, 1000*t_write/number  # ms
