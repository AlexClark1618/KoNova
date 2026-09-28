"""KN3 detector configuration."""

DETECTOR_NAME = 'KN3'

# Layer -> channel-number offset.
# Layers 1,2,3 ASIC are inverted, so its channel numbers range from 512-575.
OFFSETS = {1: 0, 2: 256, 3: 512, 4: 832}

# Layer -> bar/channel mapping file (in Bar_CH_Mappings/)
LAYER_MAPS = {
    1: 'bar_ch_map_inverted.csv',
    2: 'bar_ch_map_inverted.csv',
    3: 'bar_ch_map_inverted.csv',
    4: 'bar_ch_map_norm.csv',
}
