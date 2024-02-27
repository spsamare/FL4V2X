%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shifting VUE to correct lanes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [...
    TX_location_shifted, ...
    RX_location_shifted ...
    ] = Set_VUE_Lanes(...
    TX_location_current, ...
    RX_location_current, ...
    VUE_NS_direction_TX, ...
    VUE_NS_direction_RX, ...
    VUE_EW_direction_TX, ...
    VUE_EW_direction_RX ...
    )

catalyst = 10; % how quickly VUEs move to their lane after an intersection/junction

TX_location_shifted = TX_location_current;
RX_location_shifted = RX_location_current;


TX_location_shifted(VUE_NS_direction_TX==1, 1) = ...
    TX_location_current(VUE_NS_direction_TX==1, 1) + ...
    2*min( 1, catalyst*abs( sin( pi*...
    (TX_location_current(VUE_NS_direction_TX==1, 2)-120)/120) ) );
TX_location_shifted(VUE_NS_direction_TX==-1, 1) = ...
    TX_location_current(VUE_NS_direction_TX==-1, 1) - ...
    2*min( 1, catalyst*abs( sin( pi*...
    (TX_location_current(VUE_NS_direction_TX==-1, 2)-120)/120) ) );

RX_location_shifted(VUE_NS_direction_RX==1, 1) = ...
    RX_location_current(VUE_NS_direction_RX==1, 1) + ...
    2*min( 1, catalyst*abs( sin( pi*...
    (RX_location_current(VUE_NS_direction_RX==1, 2)-120)/120) ) );
RX_location_shifted(VUE_NS_direction_RX==-1, 1) = ...
    RX_location_current(VUE_NS_direction_RX==-1, 1) - ...
    2*min( 1, catalyst*abs( sin( pi*...
    (RX_location_current(VUE_NS_direction_RX==-1, 2)-120)/120) ) );

TX_location_shifted(VUE_EW_direction_TX==1, 2) = ...
    TX_location_current(VUE_EW_direction_TX==1, 2) - ...
    2*min( 1, catalyst*abs( sin( pi*...
    (TX_location_current(VUE_EW_direction_TX==1, 1)-120)/120) ) );
TX_location_shifted(VUE_EW_direction_TX==-1, 2) = ...
    TX_location_current(VUE_EW_direction_TX==-1, 2) + ...
    2*min( 1, catalyst*abs( sin( pi*...
    (TX_location_current(VUE_EW_direction_TX==-1, 1)-120)/120) ) );

RX_location_shifted(VUE_EW_direction_RX==1, 2) = ...
    RX_location_current(VUE_EW_direction_RX==1, 2) - ...
    2*min( 1, catalyst*abs( sin( pi*...
    (TX_location_current(VUE_EW_direction_RX==1, 1)-120)/120) ) );
RX_location_shifted(VUE_EW_direction_RX==-1, 2) = ...
    RX_location_current(VUE_EW_direction_RX==-1, 2) + ...
    2*min( 1, catalyst*abs( sin( pi*...
    (TX_location_current(VUE_EW_direction_RX==-1, 1)-120)/120) ) );
end