//
// Created by Dennis Goldfarb on 8/16/16.
//

#ifndef MSFRAGMENTDEISOTOPER_TENSORSPLINEMODEL_H
#define MSFRAGMENTDEISOTOPER_TENSORSPLINEMODEL_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>

#include "Base64.h"

struct ModelAttributes {
    unsigned int num_sulfur;
    unsigned int num_comp_sulfur;
    unsigned int num_selenium;
    unsigned int num_comp_selenium;
    unsigned int precursor_isotope;
    unsigned int fragment_isotope;

    ModelAttributes() {};

    ModelAttributes(unsigned int num_sulfur, unsigned int num_comp_sulfur, unsigned int num_selenium,
                    unsigned int num_comp_selenium, unsigned int precursor_isotope, unsigned int fragment_isotope) :
            num_sulfur(num_sulfur), num_comp_sulfur(num_comp_sulfur), num_selenium(num_selenium),
            num_comp_selenium(num_comp_selenium), precursor_isotope(precursor_isotope), fragment_isotope(fragment_isotope) {};
};

bool operator<(const ModelAttributes & lhs, const ModelAttributes & rhs);

class TensorSplineModel {

private :

    unsigned int get_xml_int(xercesc::DOMNamedNodeMap* atts, std::string att_name);
    void parse_model_attributes(xercesc::DOMNamedNodeMap* atts);
    void parse_breaks_fragment_mass(xercesc::DOMNode* current_node);
    void parse_breaks_precursor_mass(xercesc::DOMNode* current_node);
    void parse_coefficients(xercesc::DOMNode* current_node);

public :
    TensorSplineModel() {};

    TensorSplineModel(xercesc::DOMNode* current_node, xercesc::DOMNodeIterator* walker) {
        parse_xml_model(current_node, walker);
    };

    ~TensorSplineModel() {

        int max_j = 0;
        for (int i = 0; i < num_breaks_precursor_mass; ++i) {
            for (; max_j < num_breaks_fragment_mass && breaks_fragment_mass[max_j] <= breaks_precursor_mass[i]; ++max_j) {}
            for (int j = 0; j < max_j; ++j) {
                delete[] coefficients[i][j];
            }
            delete[] coefficients[i];
        }
        delete[] coefficients;

        delete[] breaks_fragment_mass;
        delete[] breaks_precursor_mass;

    }

    void parse_xml_model(xercesc::DOMNode* current_node, xercesc::DOMNodeIterator* walker);

    ModelAttributes get_model_attributes();

    float evaluate_model(float precursor_mass, float fragment_mass, bool verbose);

    unsigned int num_sulfur;
    unsigned int num_comp_sulfur;
    unsigned int num_selenium;
    unsigned int num_comp_selenium;

    unsigned int precursor_isotope;
    unsigned int fragment_isotope;

    unsigned int order_fragment_mass;
    unsigned int order_precursor_mass;

    unsigned int num_breaks_fragment_mass;
    unsigned int num_breaks_precursor_mass;

    float* breaks_fragment_mass;
    float* breaks_precursor_mass;
    float*** coefficients;
};




#endif //MSFRAGMENTDEISOTOPER_TENSORSPLINEMODEL_H
