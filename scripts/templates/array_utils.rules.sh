echo " prefix='i4' type='INTEGER(int_4)' type_temp='INTEGER(int_4)' suffix='r1' rank=':'   dim=''    alloc_bounds='my_lb:my_lb+n-1'                                    copy_bounds='1:ub'             check_bounds='SIZE( this ) < n'         get_bounds='ub = MIN( n, UBOUND( this, 1 ) )'    new_bounds='nb=MAX( n, UBOUND( this, 1 ) )' "
echo " prefix='i4' type='INTEGER(int_4)' type_temp='INTEGER(int_4)' suffix='r2' rank=':,:' dim='(2)' alloc_bounds='my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1' copy_bounds='1:ub(1), 1:ub(2)' check_bounds='ANY( SHAPE( this ) < n )' get_bounds='ub(:) = MIN( n(:), UBOUND( this ) )' new_bounds='nb(:)=MAX( n(:), UBOUND( this ) )' "


echo " prefix='l'  type='LOGICAL'        type_temp='LOGICAL'        suffix='r1' rank=':'   dim=''    alloc_bounds='my_lb:my_lb+n-1'                                    copy_bounds='1:ub'             check_bounds='SIZE( this ) < n'         get_bounds='ub = MIN( n, UBOUND( this, 1 ) )'    new_bounds='nb=MAX( n, UBOUND( this, 1 ) )' "
echo " prefix='l'  type='LOGICAL'        type_temp='LOGICAL'        suffix='r2' rank=':,:' dim='(2)' alloc_bounds='my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1' copy_bounds='1:ub(1), 1:ub(2)' check_bounds='ANY( SHAPE( this ) < n )' get_bounds='ub(:) = MIN( n(:), UBOUND( this ) )' new_bounds='nb(:)=MAX( n(:), UBOUND( this ) )' "



# 'suffix'      : [ 'r1', 'r2', 'r3', 'r4'],
# 'rank'        : [ ':', ':,:', ':,:,:', ':,:,:,:'],
# 'dim'         : [ '', '(2)', '(3)', '(4)'],
# 'alloc_bounds': [ 'my_lb:my_lb+n-1', 'my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1', 'my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1, my_lb
#(3):my_lb(3)+n(3)-1', 'my_lb(1):my_lb(1)+n(1)-1, my_lb(2):my_lb(2)+n(2)-1, my_lb(3):my_lb(3)+n(3)-1, my_lb(4):my_lb(4)+n(4)-1'],
# 'copy_bounds' : [ '1:ub', '1:ub(1), 1:ub(2)', '1:ub(1), 1:ub(2), 1:ub(3)', '1:ub(1), 1:ub(2), 1:ub(3), 1:ub(4)'],
# 'check_bounds': ['SIZE( this ) < n', 'ANY( SHAPE( this ) < n )', 'ANY( SHAPE( this ) < n )', 'ANY( SHAPE( this ) < n )'],
# 'get_bounds'  : ['ub = MIN( n, UBOUND( this, 1 ) )', 'ub(:) = MIN( n(:), UBOUND( this ) )', 'ub(:) = MIN( n(:), UBOUND( this ) )', 'ub(:) = MIN( n(:), UBOUND( this ) )']


