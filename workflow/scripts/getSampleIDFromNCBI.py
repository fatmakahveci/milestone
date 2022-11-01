#!/usr/bin/env python

from selenium import webdriver
from selenium.webdriver.common.keys import Keys

if __name__ == "__main__":
    try:
        browser = webdriver.Chrome()
        browser.get("https://www.ncbi.nlm.nih.gov/sra/?term=%22Acinetobacter+baumannii%22%5Borgn%5D")
        browser.implicitly_wait(3)
        html = browser.page_source
        results = browser.find_elements_by_class_name("rslt")
        results_texts = list()

        id_file = open("A_baumannii_id.txt", "w")

        for result in results:
            result_text = result.text
            id = result_text.index("Accession:")
            id_text = result_text[id + len("Accession: "):].strip().replace('X', 'R', 2)
            results_texts.append(id_text)

        for i in range(2,371): # change page number
            page_element = browser.find_element_by_id("pageno")
            page_element.click()
            page_element.clear()
            pageno = str(i)
            page_element.send_keys(pageno)
            browser.implicitly_wait(3)
            page_element.send_keys(Keys.ENTER)
            browser.implicitly_wait(3)
            html = browser.page_source
            results = browser.find_elements_by_class_name("rslt")

            for result in results:
                result_text = result.text
                id = result_text.index("Accession:")
                id_text = result_text[id+len("Accession: "):].strip().replace('X', 'R', 2)
                results_texts.append(id_text)

        for result in results_texts:
            id_file.write(str(result)+"\n")

    finally:
        id_file.close()
        browser.quit()
