select case (nr)
case(1)
model % coeff => COEFF_SS73DYF
model % ny = 3
model % nbl = 1
model % nbr = 2
model % coeff_L => COEFF_SS73DYF_BL
model % coeff_R => COEFF_SS73DYF_BR
model % ix = [1,2,2,3,-99,-99]
case(2)
model % coeff => COEFF_SS73COR
model % ny = 4
model % nbl = 2
model % nbr = 2
model % coeff_L => COEFF_SS73COR_BL
model % coeff_R => COEFF_SS73COR_BR
model % ix = [1,2,4,3,-99,-99]
case(3)
model % coeff => COEFF_MAGNDYF
model % ny = 4
model % nbl = 2
model % nbr = 2
model % coeff_L => COEFF_MAGNDYF_BL
model % coeff_R => COEFF_MAGNDYF_BR
model % ix = [1,2,2,3,4,-99]
case(4)
model % coeff => COEFF_MAGNCOR
model % ny = 5
model % nbl = 3
model % nbr = 2
model % coeff_L => COEFF_MAGNCOR_BL
model % coeff_R => COEFF_MAGNCOR_BR
model % ix = [1,2,4,3,5,-99]
case(5)
model % coeff => COEFF_SS73DYFCND
model % ny = 4
model % nbl = 2
model % nbr = 2
model % coeff_L => COEFF_SS73DYFCND_BL
model % coeff_R => COEFF_SS73DYFCND_BR
model % ix = [1,2,2,3,-99,4]
case(6)
model % coeff => COEFF_SS73CORCND
model % ny = 5
model % nbl = 3
model % nbr = 2
model % coeff_L => COEFF_SS73CORCND_BL
model % coeff_R => COEFF_SS73CORCND_BR
model % ix = [1,2,4,3,-99,5]
case(7)
model % coeff => COEFF_MAGNDYFCND
model % ny = 5
model % nbl = 3
model % nbr = 2
model % coeff_L => COEFF_MAGNDYFCND_BL
model % coeff_R => COEFF_MAGNDYFCND_BR
model % ix = [1,2,2,3,4,5]
case(8)
model % coeff => COEFF_MAGNCORCND
model % ny = 6
model % nbl = 4
model % nbr = 2
model % coeff_L => COEFF_MAGNCORCND_BL
model % coeff_R => COEFF_MAGNCORCND_BR
model % ix = [1,2,4,3,5,6]
case default
error stop
end select
