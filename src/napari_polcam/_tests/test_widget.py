#import numpy as np

#from napari_polcam import StokesQWidget


def test_sum():
    assert sum([1, 2, 3]) == 6, "Should be 6"

def test_sum_tuple():
    assert sum((1, 2, 2)) == 6, "Should be 6"

## make_napari_viewer is a pytest fixture that returns a napari viewer object
## capsys is a pytest fixture that captures stdout and stderr output streams
#def test_example_q_widget(make_napari_viewer, ):
#    # make viewer and add an image layer using our fixture
#    viewer = make_napari_viewer()
#    viewer.add_image(np.random.random((100, 100)))
#
#    # create our widget, passing in the viewer
#    try:
#        my_widget = StokesQWidget(viewer)
#        return True
#    except:
#        return False