//
// Created by Dennis Goldfarb on 8/16/16.
//

#include "TensorSplineModel.h"

void TensorSplineModel::parse_model_attributes(xercesc::DOMNamedNodeMap *atts) {
    num_sulfur = get_xml_int(atts,"S");
    num_comp_sulfur = get_xml_int(atts,"CS");
    num_selenium = get_xml_int(atts,"Se");
    num_comp_selenium = get_xml_int(atts,"CSe");

    precursor_isotope = get_xml_int(atts,"PrecursorIndex");
    fragment_isotope = get_xml_int(atts,"FragmentIndex");
    order_fragment_mass = get_xml_int(atts,"OrderX");
    order_precursor_mass = get_xml_int(atts,"OrderY");
}

unsigned int TensorSplineModel::get_xml_int(xercesc::DOMNamedNodeMap* atts, std::string att_name) {
    ::XMLCh* att = xercesc::XMLString::transcode(att_name.c_str());
    unsigned int ret = xercesc::XMLString::parseInt(atts->getNamedItem(att)->getNodeValue());
    xercesc::XMLString::release(&att);
    return ret;
}

void TensorSplineModel::parse_breaks_fragment_mass(xercesc::DOMNode* current_node) {
    char* input_c = xercesc::XMLString::transcode(current_node->getTextContent());
    std::string input(input_c);
    int inputLength = input.length();
    xercesc::XMLString::release(&input_c);

    int decodedLength = base64_dec_len(&input[0], inputLength);
    char decoded[decodedLength];
    base64_decode(decoded, &input[0], inputLength);

    xercesc::DOMNamedNodeMap* atts = current_node->getAttributes();
    ::XMLCh* length_att = xercesc::XMLString::transcode("length");
    int numElements = xercesc::XMLString::parseInt(atts->getNamedItem(length_att)->getNodeValue());
    xercesc::XMLString::release(&length_att);

    num_breaks_fragment_mass = numElements;
    breaks_fragment_mass = new float[numElements];

    for (int i = 0; i < numElements; i++) {
        breaks_fragment_mass[i] = *(float*) (decoded+(i*sizeof(float)));
    }
}

void TensorSplineModel::parse_breaks_precursor_mass(xercesc::DOMNode* current_node) {
    char* input_c = xercesc::XMLString::transcode(current_node->getTextContent());
    std::string input(input_c);
    int inputLength = input.length();
    xercesc::XMLString::release(&input_c);

    int decodedLength = base64_dec_len(&input[0], inputLength);
    char decoded[decodedLength];
    base64_decode(decoded, &input[0], inputLength);

    xercesc::DOMNamedNodeMap* atts = current_node->getAttributes();
    ::XMLCh* length_att = xercesc::XMLString::transcode("length");
    int numElements = xercesc::XMLString::parseInt(atts->getNamedItem(length_att)->getNodeValue());
    xercesc::XMLString::release(&length_att);
    num_breaks_precursor_mass = numElements;
    breaks_precursor_mass = new float[numElements];

    for (int i = 0; i < numElements; i++) {
        breaks_precursor_mass[i] = *(float*) (decoded+(i*sizeof(float)));
    }
}

void TensorSplineModel::parse_coefficients(xercesc::DOMNode* current_node) {
    char* input_c = xercesc::XMLString::transcode(current_node->getTextContent());
    std::string input(input_c);
    int inputLength = input.length();
    xercesc::XMLString::release(&input_c);

    int decodedLength = base64_dec_len(&input[0], inputLength);
    char decoded[decodedLength];
    base64_decode(decoded, &input[0], inputLength);


    xercesc::DOMNamedNodeMap* atts = current_node->getAttributes();
    ::XMLCh* length_att = xercesc::XMLString::transcode("length");
    int numElements = xercesc::XMLString::parseInt(atts->getNamedItem(length_att)->getNodeValue());
    xercesc::XMLString::release(&length_att);

    int num_coefficients = order_fragment_mass * order_precursor_mass;

    // allocate memory
    coefficients = new float**[num_breaks_precursor_mass];
    int max_j = 0;
    for (int i = 0; i < num_breaks_precursor_mass; ++i) {
        for (; max_j < num_breaks_fragment_mass && breaks_fragment_mass[max_j] < breaks_precursor_mass[i]; ++max_j) {}
        coefficients[i] = new float*[max_j+1];
        for (int j = 0; j <= max_j; ++j) {
            coefficients[i][j] = new float[num_coefficients];
        }
    }


    int p_index=0, f_index=-1, c_index=0, last_f_index=0;
    float precursor_mass = breaks_precursor_mass[p_index], fragment_mass;

    /*std::cout << "--MODEL DESCRIPTION--" << std::endl;
    std::cout << "S: " << num_sulfur << " CS: " << num_comp_sulfur << " Precursor isotope: "
    << precursor_isotope << " Fragment isotope: " << fragment_isotope << std::endl;*/

    for (int i = 0; i < numElements; i++) {
        c_index = i%num_coefficients;

        if (c_index==0) {
            f_index++;
            fragment_mass = breaks_fragment_mass[f_index];

            if (fragment_mass > precursor_mass || f_index == num_breaks_fragment_mass) {
                f_index=0;
                p_index++;
                precursor_mass = breaks_precursor_mass[p_index];
            }

            /*std::cout << std::endl << "Precursor mass: " << breaks_precursor_mass[p_index]
            << " Fragment mass: " << fragment_mass << " Precursor index: " << p_index
            << " Fragment index: " << f_index << std::endl;*/
        }

        coefficients[p_index][f_index][c_index] = *(float*) (decoded+(i*sizeof(float)));
        /*if (std::isnan(coefficients[p_index][f_index][c_index]) && precursor_isotope == 1) {
            std::cout << "--MODEL DESCRIPTION--" << std::endl;
            std::cout << "S: " << num_sulfur << " CS: " << num_comp_sulfur << " Precursor isotope: "
                        << precursor_isotope << " Fragment isotope: " << fragment_isotope << std::endl;
            std::cout << p_index << " " << f_index << " " << c_index << std::endl;
        }*/

        /*std::cout << "C" << c_index << ": " << coefficients[p_index][f_index][c_index] << "  \t";*/
    }

}

void TensorSplineModel::parse_xml_model(xercesc::DOMNode *current_node, xercesc::DOMNodeIterator* walker) {

    parse_model_attributes(current_node->getAttributes());


    for (current_node = walker->nextNode(); current_node != 0; current_node = walker->nextNode()) {

        char* thisNodeName = xercesc::XMLString::transcode(current_node->getNodeName());
        std::string name(thisNodeName);

        if (name == "breaksFragmentMassArrayBinary") {
            parse_breaks_fragment_mass(current_node);
        } else if (name == "breaksPrecursorMassArrayBinary") {
            parse_breaks_precursor_mass(current_node);
        } else if (name == "coefficientsArrayBinary") {
            parse_coefficients(current_node);
            xercesc::XMLString::release(&thisNodeName);
            return;
        }

        xercesc::XMLString::release(&thisNodeName);
    }
}

float TensorSplineModel::evaluate_model(float pmass, float fmass, bool verbose) {
    if (fmass >= pmass) return -1;
    if (fmass <= breaks_fragment_mass[0] || fmass >= breaks_fragment_mass[num_breaks_fragment_mass-2]) return -1;
    if (pmass <= breaks_precursor_mass[0] || pmass >= breaks_precursor_mass[num_breaks_precursor_mass-2]) return -1;

    // find index in precursor breaks
    int precursor_index = std::lower_bound(breaks_precursor_mass, breaks_precursor_mass+num_breaks_precursor_mass, pmass)-breaks_precursor_mass-1;

    // find index in fragment breaks
    int fragment_index = std::lower_bound(breaks_fragment_mass, breaks_fragment_mass+num_breaks_fragment_mass, fmass)-breaks_fragment_mass-1;

    //if (breaks_fragment_mass[fragment_index] >= breaks_precursor_mass[precursor_index+1]) return -1;

    // do the math
    float* c = coefficients[precursor_index][fragment_index];

    float x = fmass-breaks_fragment_mass[fragment_index];
    float x2 = x*x;
    float x3 = x2*x;

    float y = pmass-breaks_precursor_mass[precursor_index];
    float y2 = y*y;
    float y3 = y2*y;

    float v1 = c[15] + c[14]*y + c[13]*y2 + c[12]*y3;
    float v2 = c[11]*x + c[10]*x*y + c[9]*x*y2 + c[8]*x*y3;
    float v3 = c[7]*x2 + c[6]*x2*y + c[5]*x2*y2 + c[4]*x2*y3;
    float v4 = c[3]*x3 + c[2]*x3*y + c[1]*x3*y2 + c[0]*x3*y3;

    float v = v1+v2+v3+v4;

    if (verbose) {
        std::cout << "S: " << num_sulfur << " CS: " << num_comp_sulfur << " Precursor isotope: "
        << precursor_isotope << " Fragment isotope: " << fragment_isotope << std::endl << std::endl;

        std::cout << "Precursor index: " << precursor_index << " Fragment index: "
        << fragment_index << " precursor mass: " << pmass << " " << " fragment mass: " << fmass << std::endl << std::endl;

        std::cout << "x: " << x << " y: " << y <<  " v1: " << v1 << " v2: " << v2
        << " v3: " << v3 << " v4: " << v4 << " v: " << v << std::endl << std::endl;

        for (int i = 0; i < 16; i++) {
            std::cout << "C" << i << ": " << c[i] << "  \t";
            if (i%4 == 3) std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    return v;
}

bool operator<(const ModelAttributes &lhs, const ModelAttributes &rhs) {
    if (lhs.precursor_isotope != rhs.precursor_isotope) return lhs.precursor_isotope < rhs.precursor_isotope;
    if (lhs.fragment_isotope != rhs.fragment_isotope) return lhs.fragment_isotope < rhs.fragment_isotope;
    if (lhs.num_sulfur != rhs.num_sulfur) return lhs.num_sulfur < rhs.num_sulfur;
    if (lhs.num_comp_sulfur != rhs.num_comp_sulfur) return lhs.num_comp_sulfur < rhs.num_comp_sulfur;
    if (lhs.num_selenium != rhs.num_selenium) return lhs.num_selenium < rhs.num_selenium;
    return lhs.num_comp_selenium < rhs.num_comp_selenium;
}

ModelAttributes TensorSplineModel::get_model_attributes() {
    return ModelAttributes(num_sulfur, num_comp_sulfur, num_selenium, num_comp_selenium, precursor_isotope, fragment_isotope);
}
