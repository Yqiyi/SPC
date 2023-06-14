#include "PolarCode.h"
#include <iostream>
#include <cmath>        /* log */
#include <sstream>      // std::stringstream
#include <fstream>
#include <iomanip>      // std::setprecision
#include <random>
#include <algorithm>
#include <chrono>
#include <vector>
#include <functional>
using namespace std;


void PolarCode::initialize_frozen_bits() {
    std::vector<double> channel_vec(_block_length);
	uint64_t msg_i = 0;
    for (uint64_t i = 0; i < _block_length; ++i) {
        channel_vec.at(i) = _design_epsilon;
    }
    for (uint8_t iteration = 0; iteration < _n; ++iteration) {
        uint64_t  increment = 1 << iteration;
        for (uint64_t j = 0; j < increment; j +=  1) {
            for (uint64_t i = 0; i < _block_length; i += 2 * increment) {
                double c1 = channel_vec.at(i + j);
                double c2 = channel_vec.at(i + j + increment);
                channel_vec.at(i + j) = c1 + c2 - c1 * c2;
                channel_vec.at(i + j + increment) = c1 * c2;
            }
        }
    }

    _channel_order_descending.resize(_block_length);
    std::size_t n_t(0);
    std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
    std::sort(  std::begin(_channel_order_descending),
                std::end(_channel_order_descending),
                [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] < channel_vec[_bit_rev_order.at(i2)]; } );

    uint64_t  effective_info_length = _info_length + _crc_size;

    for (uint64_t i = 0; i < effective_info_length; ++i) {
        _frozen_bits.at(_channel_order_descending.at(i)) = 2;       //information bit = 2
    }
    for (uint64_t i = effective_info_length; i < _block_length; ++i) {
        _frozen_bits.at(_channel_order_descending.at(i)) = _msg_bits.at(msg_i++);      //frozen bit = 0/1 as msg_bits
    }
}

// extract message by polar coding
std::vector<uint8_t> PolarCode::encode_msg_ext(std::vector<uint8_t> info_bits_padded) {

	std::vector<uint8_t> coded_bits(_block_length);
	std::vector<uint8_t> msg_bits_ext(_msg_length);
	uint64_t msg_i = 0;

	for (uint8_t iteration = 0; iteration < _n; ++iteration) {
		uint64_t  increment = (uint64_t)(1 << iteration);
		for (uint64_t j = 0; j < increment; j += 1) {
			for (uint64_t i = 0; i < _block_length; i += 2 * increment) {
				info_bits_padded.at(i + j) = (uint8_t)(info_bits_padded.at(i + j) + info_bits_padded.at(i + j + increment));
			}
		}
	}

	for (uint64_t i = 0; i < _block_length; ++i) {
		coded_bits.at(i) = (info_bits_padded.at(_bit_rev_order.at(i)) %= 2);
	}
	// extract message with frozen indices
	for (uint64_t i = _info_length; i < _block_length; ++i) {
		msg_bits_ext.at(msg_i++) = coded_bits.at(_channel_order_descending.at((i)));      // frozen bit = 0/1 as msg_bits
	}
	return msg_bits_ext;
}

// generate stego by polar coding
std::vector<uint8_t> PolarCode::encode_msg(std::vector<uint8_t> info_bits_padded) {

	std::vector<uint8_t> coded_bits(_block_length);
	for (uint8_t iteration = 0; iteration < _n; ++iteration) {
		uint64_t  increment = (uint64_t)(1 << iteration);
		for (uint64_t j = 0; j < increment; j += 1) {
			for (uint64_t i = 0; i < _block_length; i += 2 * increment) {
				info_bits_padded.at(i + j) = (uint8_t)(info_bits_padded.at(i + j) + info_bits_padded.at(i + j + increment));
			}
		}
	}

	for (uint64_t i = 0; i < _block_length; ++i) {
		coded_bits.at(i) = (info_bits_padded.at(_bit_rev_order.at(i)) %= 2);
	}

	return coded_bits;
}

bool PolarCode::crc_check(uint8_t * info_bit_padded) {
    bool crc_pass = true;
    for (uint64_t i = _info_length; i < _info_length + _crc_size; ++i) {
        uint8_t  crc_bit = 0;
        for (uint64_t j = 0; j < _info_length; ++j) {
            crc_bit = (uint8_t) ((crc_bit + _crc_matrix.at(i - _info_length).at(j) * info_bit_padded[_channel_order_descending.at(j)]) % 2);
        }

        if (crc_bit != info_bit_padded[_channel_order_descending.at(i)]) {
            crc_pass = false;
            break;
        }
    }

    return crc_pass;
}

std::vector<uint8_t> PolarCode::decode_scl_p1(std::vector<double> p1, std::vector<double> p0, uint64_t list_size) {

    _list_size = list_size;
    _llr_based_computation = false;

    initializeDataStructures();

    uint64_t  l = assignInitialPath();

    double * p_0 = getArrayPointer_P(0, l);

    for (uint64_t beta = 0; beta < _block_length; ++beta ) {
        p_0[2*beta] = (double) p0.at(beta);
        p_0[2*beta + 1] = (double) p1.at(beta);
    }

    return decode_scl();

}

std::vector<uint8_t> PolarCode::decode_scl_llr(std::vector<double> llr, uint64_t list_size) {

    _list_size = list_size;

    _llr_based_computation = true;

    initializeDataStructures();

    uint64_t  l = assignInitialPath();

    double * llr_0 = getArrayPointer_LLR(0, l);

    for (uint64_t beta = 0; beta < _block_length; ++beta ) {
        llr_0[beta] = llr.at(beta);
    }

    return decode_scl();

}

std::vector<uint8_t> PolarCode::decode_scl() {

    for (uint64_t phi = 0; phi < _block_length; ++phi ){

        if (_llr_based_computation )
            recursivelyCalcLLR(_n, phi);
        else
            recursivelyCalcP(_n, phi);


        if (_frozen_bits.at(phi) != 2)      // frozen bit = 0/1 as msg_bits
            continuePaths_FrozenBit(phi);
        else
            continuePaths_UnfrozenBit(phi);

        if ((phi%2) == 1)
            recursivelyUpdateC(_n, phi);

    }
    //uint64_t l = findMostProbablePath((bool) _crc_size);
	uint64_t l = 0;

    uint8_t * c_0 = _arrayPointer_Info.at(l);
    std::vector<uint8_t> deocded_bits(_block_length);
    for (uint64_t beta = 0; beta < _block_length; ++beta )
        deocded_bits.at(beta) = c_0[beta];

    for (uint64_t s = 0; s < _list_size; ++s) {
        delete[] _arrayPointer_Info.at(s);
        for (uint64_t lambda = 0; lambda < _n + 1; ++lambda) {

            if (_llr_based_computation )
                delete[] _arrayPointer_LLR.at(lambda).at(s);
            else
                delete[] _arrayPointer_P.at(lambda).at(s);
            delete[] _arrayPointer_C.at(lambda).at(s);
        }
    }

    return deocded_bits;

}


// the following are subfunctions for SCL decoding
void PolarCode::initializeDataStructures() {

    while (_inactivePathIndices.size()) {
        _inactivePathIndices.pop();
    };
    _activePath.resize(_list_size);

    if (_llr_based_computation) {
        _pathMetric_LLR.resize(_list_size);
        _arrayPointer_LLR.resize(_n + 1);
        for (int i = 0; i < _n + 1; ++i)
            _arrayPointer_LLR.at(i).resize(_list_size);
    }
    else {
        _arrayPointer_P.resize(_n + 1);
        for (int i = 0; i < _n + 1; ++i)
            _arrayPointer_P.at(i).resize(_list_size);
    }

    _arrayPointer_C.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayPointer_C.at(i).resize(_list_size);

    _arrayPointer_Info.resize(_list_size);

    _pathIndexToArrayIndex.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _pathIndexToArrayIndex.at(i).resize(_list_size);

    _inactiveArrayIndices.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i) {
        while (_inactiveArrayIndices.at(i).size()) {
            _inactiveArrayIndices.at(i).pop();
        };
    }

    _arrayReferenceCount.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayReferenceCount.at(i).resize(_list_size);

    for (uint64_t s = 0; s < _list_size; ++s) {
        _arrayPointer_Info.at(s) = new uint8_t[_block_length]();
        for (uint64_t lambda = 0; lambda < _n + 1; ++lambda) {
            if (_llr_based_computation) {
                _arrayPointer_LLR.at(lambda).at(s) = new double[(1 << (_n - lambda))]();
            }
            else {
                _arrayPointer_P.at(lambda).at(s) = new double[2 * (1 << (_n - lambda))]();
            }
            _arrayPointer_C.at(lambda).at(s) = new uint8_t[2 * (1 << (_n - lambda))]();
            _arrayReferenceCount.at(lambda).at(s) = 0;
            _inactiveArrayIndices.at(lambda).push(s);
        }
    }

    for (uint64_t l = 0; l < _list_size; ++l) {
        _activePath.at(l) = 0;
        _inactivePathIndices.push(l);
        if (_llr_based_computation) {
            _pathMetric_LLR.at(l) = 0;
        }
    }
}

uint64_t PolarCode::assignInitialPath() {

    uint64_t  l = _inactivePathIndices.top();
    _inactivePathIndices.pop();
    _activePath.at(l) = 1;
    // Associate arrays with path index
    for (uint64_t lambda = 0; lambda < _n + 1; ++lambda) {
        uint64_t  s = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();
        _pathIndexToArrayIndex.at(lambda).at(l) = s;
        _arrayReferenceCount.at(lambda).at(s) = 1;
    }
    return l;
}

uint64_t PolarCode::clonePath(uint64_t l) {
    uint64_t l_p = _inactivePathIndices.top();
    _inactivePathIndices.pop();
    _activePath.at(l_p) = 1;

    if (_llr_based_computation)
        _pathMetric_LLR.at(l_p) = _pathMetric_LLR.at(l);

    for (uint64_t lambda = 0; lambda < _n + 1; ++lambda ) {
        uint64_t s = _pathIndexToArrayIndex.at(lambda).at(l);
        _pathIndexToArrayIndex.at(lambda).at(l_p) = s;
        _arrayReferenceCount.at(lambda).at(s)++;
    }
    return l_p;
}

void PolarCode::killPath(uint64_t l) {
    _activePath.at(l) = 0;
    _inactivePathIndices.push(l);
    if (_llr_based_computation )
        _pathMetric_LLR.at(l) = 0;

    for (uint64_t lambda = 0; lambda < _n + 1; ++lambda ) {
        uint64_t s = _pathIndexToArrayIndex.at(lambda).at(l);
        _arrayReferenceCount.at(lambda).at(s)--;
        if (_arrayReferenceCount.at(lambda).at(s) == 0 ) {
            _inactiveArrayIndices.at(lambda).push(s);
        }
    }
}

double * PolarCode::getArrayPointer_P(uint64_t lambda, uint64_t  l) {
    uint64_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint64_t s_p;
    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {
        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_P.at(lambda).at(s_p));
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;
    }
    return _arrayPointer_P.at(lambda).at(s_p);
}

double * PolarCode::getArrayPointer_LLR(uint64_t lambda, uint64_t  l) {
    uint64_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint64_t s_p;
    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {
        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));
        std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) +  (1 << (_n - lambda)),  _arrayPointer_LLR.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;
    }
    return _arrayPointer_LLR.at(lambda).at(s_p);
}


uint8_t * PolarCode::getArrayPointer_C(uint64_t lambda, uint64_t  l) {
    uint64_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint64_t s_p;
    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {

        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        if (_llr_based_computation )
            std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) +  (1 << (_n - lambda)),  _arrayPointer_LLR.at(lambda).at(s_p));
        else
            std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_P.at(lambda).at(s_p));

        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;

    }
    return _arrayPointer_C.at(lambda).at(s_p);
}

void PolarCode::recursivelyCalcP(uint64_t lambda, uint64_t phi) {
    if ( lambda == 0 )
        return;
    uint64_t psi = phi >> 1;
    if ( (phi % 2) == 0)
        recursivelyCalcP(lambda -1, psi);

    double sigma = 0.0f;
    for (uint64_t l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        double * p_lambda = getArrayPointer_P(lambda, l);
        double * p_lambda_1 = getArrayPointer_P(lambda - 1, l);

        uint8_t * c_lambda = getArrayPointer_C(lambda, l);
        for (uint64_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            if ( (phi %2) == 0 ){
                p_lambda[2 * beta] = 0.5f * ( p_lambda_1[2*(2*beta)]*p_lambda_1[2*(2*beta+1)]
                                              + p_lambda_1[2*(2*beta) + 1]*p_lambda_1[2*(2*beta+1) + 1]);
                p_lambda[2 * beta + 1] = 0.5f * ( p_lambda_1[2*(2*beta) +1]*p_lambda_1[2*(2*beta+1)]
                                                  + p_lambda_1[2*(2*beta)]*p_lambda_1[2*(2*beta+1) + 1]);
            }
            else {
                uint8_t  u_p = c_lambda[2*beta];
                p_lambda[2 * beta] = 0.5f * p_lambda_1[2*(2*beta) + (u_p % 2)] *   p_lambda_1[2*(2*beta + 1)];
                p_lambda[2 * beta + 1] = 0.5f * p_lambda_1[2*(2*beta) + ((u_p+1) % 2)] *   p_lambda_1[2*(2*beta + 1) + 1];
            }
            sigma = std::max(sigma,  p_lambda[2 * beta]);
            sigma = std::max(sigma,  p_lambda[2 * beta + 1]);


        }
    }

    for (uint64_t l = 0; l < _list_size; ++l) {
        if (sigma == 0) // Typically happens because of undeflow
            break;
        if (_activePath.at(l) == 0)
            continue;
        double *p_lambda = getArrayPointer_P(lambda, l);
        for (uint64_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            p_lambda[2 * beta] = p_lambda[2 * beta] / sigma;
            p_lambda[2 * beta + 1] = p_lambda[2 * beta + 1] / sigma;
        }
    }
}

void PolarCode::recursivelyCalcLLR(uint64_t lambda, uint64_t phi) {
    if ( lambda == 0 )
        return;
    uint64_t psi = phi >> 1;
    if ( (phi % 2) == 0)
        recursivelyCalcLLR(lambda -1, psi);

    for (uint64_t l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        double * llr_lambda = getArrayPointer_LLR(lambda, l);
        double * llr_lambda_1 = getArrayPointer_LLR(lambda - 1, l);

        uint8_t * c_lambda = getArrayPointer_C(lambda, l);
        for (uint64_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            if ( (phi %2) == 0 ){
                if (40 > std::max(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]))){
                    llr_lambda[beta] = std::log ( (exp(llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1]) + 1) /
                                                  (exp(llr_lambda_1[2*beta]) + exp(llr_lambda_1[2*beta+1])));
                }
                else {
                    llr_lambda[beta] = (double)  ((llr_lambda_1[2 * beta] < 0) ? -1 : (llr_lambda_1[2 * beta] > 0)) *
                                       ((llr_lambda_1[2 * beta + 1] < 0) ? -1 : (llr_lambda_1[2 * beta + 1] > 0)) *
                                       std::min( std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]));
                }
            }
            else {
                uint8_t  u_p = c_lambda[2*beta];
                llr_lambda[beta] = (1 - 2 * u_p) * llr_lambda_1[2*beta] + llr_lambda_1[2*beta + 1];
            }

        }
    }
}

void PolarCode::recursivelyUpdateC(uint64_t lambda, uint64_t phi) {

    uint64_t psi = phi >> 1;
    for (uint64_t l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        uint8_t *c_lambda = getArrayPointer_C(lambda, l);
        uint8_t *c_lambda_1 = getArrayPointer_C(lambda - 1, l);
        for (uint64_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            c_lambda_1[2 * (2 * beta) + (psi % 2)] = (uint8_t) ((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
            c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
        }
    }
    if ( (psi % 2) == 1)
        recursivelyUpdateC((uint64_t) (lambda - 1), psi);

}

void PolarCode::continuePaths_FrozenBit(uint64_t phi) {
    for (uint64_t l = 0; l < _list_size; ++ l) {
        if (_activePath.at(l) == 0)
            continue;
        uint8_t  * c_m = getArrayPointer_C(_n, l);
        c_m[(phi % 2)] = _frozen_bits.at(phi);           // frozen value assumed to be zero
        if (_llr_based_computation) {
            double *llr_p = getArrayPointer_LLR(_n, l);
			int dir = 2 * _frozen_bits.at(phi) - 1;
            _pathMetric_LLR.at(l) += log(1 +  exp(dir *llr_p[0]));     //for SPC
        }
        _arrayPointer_Info.at(l)[phi] = _frozen_bits.at(phi);
    }
}

void PolarCode::continuePaths_UnfrozenBit(uint64_t phi) {

    std::vector<double>  probForks((unsigned long) (2 * _list_size));
    std::vector<double> probabilities;
    std::vector<uint8_t>  contForks((unsigned long) (2 * _list_size));


    uint64_t  i = 0;
    for (unsigned l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0) {
            probForks.at(2 * l) = NAN;
            probForks.at(2 * l + 1) = NAN;
        }
        else {
            if (_llr_based_computation ) {
                double *llr_p = getArrayPointer_LLR(_n, l);
                probForks.at(2 * l) =  - (_pathMetric_LLR.at(l) + log(1 + exp(-llr_p[0])));
                probForks.at(2 * l + 1) = -  (_pathMetric_LLR.at(l) + log(1 + exp(llr_p[0])));
            }
            else {
                double *p_m = getArrayPointer_P(_n, l);
                probForks.at(2 * l) = p_m[0];
                probForks.at(2 * l + 1) = p_m[1];
            }

            probabilities.push_back(probForks.at(2 * l));
            probabilities.push_back(probForks.at(2 * l +1));

            i++;
        }
    }

    uint64_t  rho = _list_size;
    if ( (2*i) < _list_size)
        rho = (uint64_t) 2 * i;

    for (uint8_t l = 0; l < 2 * _list_size; ++l) {
        contForks.at(l) = 0;
    }
    std::sort(probabilities.begin(), probabilities.end(), std::greater<double>());

    double threshold = probabilities.at((unsigned long) (rho - 1));
    uint64_t num_paths_continued = 0;

    for (uint8_t l = 0; l < 2 * _list_size; ++l) {
        if (probForks.at(l) > threshold) {
            contForks.at(l) = 1;
            num_paths_continued++;
        }
        if (num_paths_continued == rho) {
            break;
        }
    }

    if  ( num_paths_continued < rho ) {
        for (uint8_t l = 0; l < 2 * _list_size; ++l) {
            if (probForks.at(l) == threshold) {
                contForks.at(l) = 1;
                num_paths_continued++;
            }
            if (num_paths_continued == rho) {
                break;
            }
        }
    }

    for (unsigned l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        if ( contForks.at(2 * l)== 0 && contForks.at(2 * l + 1) == 0 )
            killPath(l);
    }

    for (unsigned l = 0; l < _list_size; ++l) {
        if ( contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0 )
            continue;
        uint8_t * c_m = getArrayPointer_C(_n, l);

        if ( contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1 ) {

            c_m[(phi%2)] = 0;
            uint64_t l_p = clonePath(l);
            c_m = getArrayPointer_C(_n, l_p);
            c_m[(phi%2)] = 1;

            std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) +  phi,  _arrayPointer_Info.at(l_p));
            _arrayPointer_Info.at(l)[phi] = 0;
            _arrayPointer_Info.at(l_p)[phi] = 1;

            if (_llr_based_computation ) {
                double *llr_p = getArrayPointer_LLR(_n, l);
                _pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
                llr_p = getArrayPointer_LLR(_n, l_p);
                _pathMetric_LLR.at(l_p) += log(1 + exp(llr_p[0]));
            }

        }
        else {
            if ( contForks.at(2 * l) == 1) {
                c_m[(phi%2)] = 0;
                _arrayPointer_Info.at(l)[phi] = 0;

                if (_llr_based_computation ) {
                    double *llr_p = getArrayPointer_LLR(_n, l);
                    _pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
                }
            }
            else {
                c_m[(phi%2)] = 1;
                _arrayPointer_Info.at(l)[phi] = 1;
                if (_llr_based_computation ) {
                    double *llr_p = getArrayPointer_LLR(_n, l);
                    _pathMetric_LLR.at(l) += log(1 + exp(llr_p[0]));
                }
            }
        }
    }

}

uint64_t PolarCode::findMostProbablePath(bool check_crc) {

    uint64_t  l_p = 0;
    double p_p1 = 0;
    double p_llr = std::numeric_limits<double>::max();
    bool path_with_crc_pass = false;
    for (uint64_t l = 0; l < _list_size; ++l) {

        if (_activePath.at(l) == 0)
            continue;

        if ( (check_crc) && (! crc_check(_arrayPointer_Info.at(l))))
            continue;

        path_with_crc_pass = true;

        if (_llr_based_computation) {
            if (_pathMetric_LLR.at(l) < p_llr ) {
                p_llr = _pathMetric_LLR.at(l);
                l_p  = l;
            }
        }
        else {
            uint8_t * c_m = getArrayPointer_C(_n, l);
            double * p_m = getArrayPointer_P(_n, l);
            if ( p_p1 < p_m[c_m[1]]) {
                l_p = l;
                p_p1 = p_m[c_m[1]];
            }
        }
    }
    if ( path_with_crc_pass)
        return l_p;
    else
        return findMostProbablePath(false);
}

void PolarCode::create_bit_rev_order() {
    for (uint64_t i = 0; i < _block_length; ++i) {
        uint64_t to_be_reversed = i;
        _bit_rev_order.at(i) = (uint64_t) ((to_be_reversed & 1) << (_n - 1));
        for (uint8_t j = (uint8_t) (_n - 1); j; --j) {
            to_be_reversed >>= 1;
            _bit_rev_order.at(i) += (to_be_reversed & 1) << (j - 1);
        }
    }
}
