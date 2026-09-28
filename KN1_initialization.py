"""KN1 detector configuration."""

DETECTOR_NAME = 'KN1'

# Layer -> channel-number offset, e.g. layer 1 ch 64-127, layer 2 ch 320-383, etc.
OFFSETS = {1: 64, 2: 320, 3: 576, 4: 832}

# Layer -> bar/channel mapping file (in Bar_CH_Mappings/)
LAYER_MAPS = {
    1: 'bar_ch_map_norm.csv',
    2: 'bar_ch_map_norm.csv',
    3: 'bar_ch_map_norm.csv',
    4: 'bar_ch_map_norm.csv',
}
