# tests/test_widget.py

import pytest
from napari_polcam import StokesEstimation


# Using the napari viewer fixture
@pytest.mark.qt
def test_stokes_initialization(make_napari_viewer):
    # Create a real napari viewer instance
    viewer = make_napari_viewer()

    # Initialize the StokesEstimation widget with the real viewer
    widget = StokesEstimation(napari_viewer=viewer)

    # Check that the widget was created successfully
    assert widget is not None

