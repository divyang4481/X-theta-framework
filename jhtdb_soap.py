import os
import requests
from zeep import Client
from zeep.transports import Transport
from zeep.cache import InMemoryCache

WSDL = "http://turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL"

def make_client(timeout_s: int = 60) -> Client:
    # Proxy-friendly: respects HTTPS_PROXY / HTTP_PROXY env vars automatically.
    sess = requests.Session()
    sess.trust_env = True  # important behind corporate proxies

    transport = Transport(
        session=sess,
        timeout=timeout_s,
        cache=InMemoryCache(),
    )
    return Client(wsdl=WSDL, transport=transport)

def list_methods(client: Client):
    svc = list(client.wsdl.services.values())[0]
    port = list(svc.ports.values())[0]
    ops = sorted(port.binding._operations.keys())
    return ops

if __name__ == "__main__":
    client = make_client()
    print("Available SOAP methods:")
    for m in list_methods(client):
        print(" -", m)
