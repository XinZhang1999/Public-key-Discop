import os
from typing import Optional
import torch
# from scipy.io.wavfile import read, write
# from PIL import Image

from config import Settings, text_default_settings, image_default_settings, audio_default_settings
from model import get_model, get_feature_extractor, get_tokenizer
from utils import SingleExampleOutput, check_dir
import PRNEncryption
from PRNEncryption import PRNEncryption



def run_single_example(message, settings: Settings = text_default_settings, context: Optional[str] = None):
    E = PRNEncryption()  # Pseudorandom Public-key ECC encryption scheme (satisfied IND$-CPA property) based on admissible encoding
    S = ''
    ciphertext = E.encrypt(message)
    binary_string = ''.join(format(byte, '08b') for byte in ciphertext)
    for bit in binary_string:
        S += '0' if bit == '0' else '1'
    print("ciphertext = ", S)
    bitlength = len(S)
    print('bitlength=', bitlength)
    settings.length = 2 * bitlength
    
    if settings.algo == 'sample':
        from random_sample_cy import encode_text
    elif settings.algo in ['Discop', 'Discop_baseline']:
        from stega_cy import encode_text, decode_text
    else:
        raise NotImplementedError("`Settings.algo` must belong to {'Discop', 'Discop_baseline', 'sample'}!")

    if context is None:
        context = 'We were both young when I first saw you, I close my eyes and the flashback starts.'
    print("context = ",context)
    
    model = get_model(settings)
    tokenizer = get_tokenizer(settings)
    
    single_example_output: SingleExampleOutput = encode_text(model, tokenizer, S, context, settings)
    print(single_example_output)
    print("---------------------------------------------------------")
    if settings.algo != 'sample':
        message_encoded = S[:single_example_output.n_bits]
        message_decoded = decode_text(model, tokenizer, single_example_output.generated_ids, context, settings)
        message_decoded = message_decoded[:bitlength]
        print(message_encoded)
        print(message_decoded)
        
        integer_value = int(message_decoded, 2)
        byte_value = integer_value.to_bytes((len(message_decoded) + 7) // 8, byteorder='big')
        print(byte_value)
        m = E.decrypt(byte_value)
        print("m =",m)

if __name__ == '__main__':
    message = "Attack at 9:00"
    # Text Generation
    settings = text_default_settings
    # Can change to Llama2 by 
    # settings.model_name = llama2
    settings.model_name = 'gpt2'  
    
    settings.device = torch.device('cuda:0')
    # settings.algo = 'Discop_baseline'
    settings.algo = 'Discop'
    context = """Years later, he would find himself"""
    run_single_example(message, settings, context)