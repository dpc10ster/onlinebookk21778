TrapezoidalArea <- function( noise, signal )
{
	area = 0
	for( ns in 1 : length( signal ) ) {
		a = noise[ noise < signal[ ns ] ]
		b = noise[ noise == signal[ ns ] ]
		area = area + length( a )
		area = area + length( b ) * 0.5
	}
	area = area / length( noise ) / length( signal )
	return( area )
}


ROC_Area <- function( continuous_array, binary_array ) {
	if( length( binary_array[ binary_array == 0 ] ) == 0 ) {
		return( -999 )
	}
	normal_array <- continuous_array[ binary_array == 0 ]
	if( length( binary_array[ binary_array == 1 ] ) == 0 ) {
		return( -999 )
	}
	abnormal_array <- continuous_array[ binary_array == 1 ]

	return( TrapezoidalArea( normal_array, abnormal_array ) )
}