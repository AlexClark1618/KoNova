"""KN2 detector configuration."""

DETECTOR_NAME = 'KN2'

# Layer -> channel-number offset.
# Layer 3 ASIC is inverted, so its channel numbers range from 512-575.
OFFSETS = {1: 64, 2: 320, 3: 512, 4: 832}

# Layer -> bar/channel mapping file (in Bar_CH_Mappings/)
LAYER_MAPS = {
    1: 'bar_ch_map_norm.csv',
    2: 'bar_ch_map_norm.csv',
    3: 'bar_ch_map_inverted.csv',
    4: 'bar_ch_map_norm.csv',
}
